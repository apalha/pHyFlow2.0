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
# Define the panel geometries
    
# Cylinder (Radius =1, nPanels = 100)    
R       = 1.0   # Radius of cylinder
nPanel  = 100   # Number of panels
dPanel  = np.spacing(100) # Spacing between panel and colloc. point
theta   = np.linspace(np.pi,-np.pi,nPanel+1) # Panel polar angles, to define them.
dtheta  = theta[1]-theta[0] # Angle spacing
r       = (R + dPanel) / np.cos(dtheta/2.0) # Radial location of the panel end points

# Make the cylinder and append the parameters to a dictionary.
cylinderData = {'xPanel' : r*np.cos(theta - dtheta/2),
                'yPanel' : r*np.sin(theta - dtheta/2),
                'cmGlobal'   : np.array([0.,0.]),
                'thetaLocal' : 0.,
                'dPanel' : np.spacing(100)}
                
# For now, only a single cylinder
geometries = {'cylinder':cylinderData}

# Initalize panelBody 
panels = pHyFlow.panel.Panels(geometries=geometries)

#------------------------------------------------------------------------------                                  
                                   
                                   
#------------------------------------------------------------------------------
# Initialize the vortex-panel class
                       
startTime = time.time()                     
vortexPanel = pHyFlow.vortexPanel.VortexPanel(blobs,panels,
                                              couplingParams={'panelStrengthUpdate':'constant'})
print "\nTime to initialize the coupled vortex-panel problem: %g seconds." % (time.time() - startTime)
                         
                         
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
# Evaluate the velocity field - due to panels
                         
x,y  = np.meshgrid(np.linspace(-4,4,20),np.linspace(-2,2,20))
vortexPanel.panels.solve(vInf[0],vInf[1])
vx,vy = vortexPanel.panels.evaluateVelocity(x.flatten(),y.flatten())

#------------------------------------------------------------------------------                                   

#------------------------------------------------------------------------------
# Plot the convection

py.figure(1)
py.title('Convection of zero-valued blobs')
#py.clf()
py.quiver(x.flat,y.flat,vx+vInf[0],vy+vInf[1],scale=100)
py.scatter(blobs.x,blobs.y,c='g',s=30,edgecolor='none',label='blobs')
[py.plot(x,y,k+'-',label=geoName) for x,y,geoName,k in zip(vortexPanel.panels.xyPanelGlobal[0],
                                                            vortexPanel.panels.xyPanelGlobal[1],
                                                            vortexPanel.panels.geometryKeys,['b','g','k'])]
py.grid(True)
py.axis('scaled')
py.axis([-4,4,-2,2])
    
for i in xrange(1000):  
    if i % 50 == 0:     
        py.scatter(blobs.x,blobs.y,c='g',s=30,edgecolor='none',label='blobs')
        py.draw()
    print i
    vortexPanel.evolve()

#------------------------------------------------------------------------------