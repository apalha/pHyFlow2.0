"""
pHyFlow test of the coupled vortexPanel class

Description
-----------

Simple particle convection test


:First added:   2014-02-21  
:Last updated:  2014-02-21    
:Copyright:     Copyright (C) 2014 Lento Manickathan, **pHyFlow**
:Licence:       GNU GPL version 3 or any later version
"""

import pHyFlow
import numpy as np
import time
import pylab as py

py.ion()

#------------------------------------------------------------------------------
# Global variables

vInf = np.array([1.0,0.0]) # the free stream velocity


nTimeSteps = 20 # evolve the blobs for nTimeSteps

#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
# Initialize vortex blobs

# Define blobs
nBlobs = 4
overlap = 1.0 # the overlap ration between the blobs
h = 0.01 # the cell size to which each blob is associated to

deltaTc = 0.01 # the size of the time step
nu = 0.0

blobControlParams = {'stepRedistribution':0,'stepPopulationControl':0,\
                     'gThresholdLocal':1e-8,'gThresholdGlobal':1e-8}

# generate the vortex blobs
x,y = -np.ones(nBlobs)*5.0, np.linspace(-1.0,1.0,nBlobs)
g = np.zeros(nBlobs)

wField = (x,y,g) # the vorticity field made up of the blobs coordinates and circulation

# Generate the blobs
blobs = pHyFlow.vortex.VortexBlobs(wField,vInf,nu,deltaTc,h,overlap,
                                   blobControlParams=blobControlParams)
#------------------------------------------------------------------------------
                                   
#------------------------------------------------------------------------------
# Initilize the panels

# Define the panel data
panelData = dict(xPanel = [np.array([-1.0,-1.0,1.0,1.0])],
                 yPanel = [np.array([-1.0,1.0,1.0,-1.0])],
                 cmGlobal = [np.array([0.,0.])],
                 thetaLocal = [0.],
                 dPanel = [np.spacing(100)])

# Generete the panels
panels = pHyFlow.panel.Panels(panelData)

xPanel, yPanel = panels.xPanelGlobalCat, panels.yPanelGlobalCat
xCP,yCP = panels.xCPGlobalCat, panels.yCPGlobalCat
normX,normY = panels._Panels__norm[0], panels._Panels__norm[1]

#------------------------------------------------------------------------------                                   
                                   
                                   
#------------------------------------------------------------------------------
# Initialize the vortex-panel class
                                   
vortexPanel = pHyFlow.vortexPanel.VortexPanel(blobs,panels)

#------------------------------------------------------------------------------                                   

#------------------------------------------------------------------------------
# Plot 

py.figure()
py.scatter(blobs.x,blobs.y,c='g',edgecolor='none',label='blobs')
py.plot(xPanel, yPanel,'b.-',label='panel')
py.plot(xCP,yCP,'r.',label='colloc. points')
py.quiver(xCP,yCP,normX,normY)
py.grid()
py.axis('scaled')

#------------------------------------------------------------------------------