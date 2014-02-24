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


#nTimeSteps = 20 # evolve the blobs for nTimeSteps

#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
# Initialize vortex blobs

# Define blobs
nBlobs = 10
overlap = 1.0 # the overlap ration between the blobs
h = 0.01 # the cell size to which each blob is associated to

deltaTc = 0.01 # the size of the time step
nu = 0.0

blobControlParams = {'stepRedistribution':0,'stepPopulationControl':0,\
                     'gThresholdLocal':1e-8,'gThresholdGlobal':1e-8}

# generate the vortex blobs
x,y = -np.ones(nBlobs)*4.0, np.linspace(-1.5,1.5,nBlobs)
g = np.zeros(nBlobs)

wField = (x,y,g) # the vorticity field made up of the blobs coordinates and circulation

# Generate the blobs
blobs = pHyFlow.vortex.VortexBlobs(wField,vInf,nu,deltaTc,h,overlap,
                                   blobControlParams=blobControlParams)
#------------------------------------------------------------------------------
                                   
#------------------------------------------------------------------------------
# Initilize the panels

R = 1
nPanels = 100
dPanel  = np.spacing(100)
theta   = np.linspace(np.pi,-np.pi,nPanels+1) # Panel polar angles
dtheta  = theta[1]-theta[0] # Angle spacing
r       = (R + dPanel) / np.cos(dtheta/2.0) # Radial location of the panel end points

# Panel location

# Define the panel data
panelData = dict(xPanel = [r*np.cos(theta[:-1] - dtheta/2)],
                 yPanel = [r*np.sin(theta[:-1] - dtheta/2)],
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
                                   
vortexPanel = pHyFlow.vortexPanel.VortexPanel(blobs,panels,
                                              couplingParams={'panelStrengthUpdate':'constant'})

x,y  = np.meshgrid(np.linspace(-4,4,10),np.linspace(-2,2,10))
vortexPanel.panels.solve(vInf[0],vInf[1])
vx,vy = vortexPanel.panels.evaluateVelocity(x.flatten(),y.flatten())

#------------------------------------------------------------------------------                                   

#------------------------------------------------------------------------------
# Plot 

py.figure(1)
#py.clf()
py.quiver(x.flat,y.flat,vx+vInf[0],vy+vInf[1],scale=100)
py.scatter(blobs.x,blobs.y,c='g',edgecolor='none',label='blobs')
py.plot(xPanel, yPanel,'b.-',label='panel')
py.plot(xCP,yCP,'r.',label='colloc. points')
py.grid(True)
py.axis('scaled')
py.axis([-4,4,-2,2])
    
for i in xrange(1000):  
    if i % 50 == 0:     
        py.scatter(blobs.x,blobs.y,c='g',edgecolor='none',label='blobs')
        py.draw()
    print i
    vortexPanel.evolve()

#------------------------------------------------------------------------------