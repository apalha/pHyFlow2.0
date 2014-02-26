"""
Test Multi-Body Problem. 3 Cylinders of different size, with different
number of panels.

:First Added:   2013-11-19
"""

import pHyFlow

import numpy as np
import pylab as py
import time
py.ion()

#------------------------------------------------------------------------------
# Global variables

vInf = np.array([1.0,0.0]) # the free stream velocity

#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Initialize vortex blobs

# Define blobs
nBlobs = 20
overlap = 1.0 # the overlap ration between the blobs
h = 0.01 # the cell size to which each blob is associated to

deltaTc = 0.01 # the size of the time step
nu = 0.0

blobControlParams = {'stepRedistribution':0,'stepPopulationControl':0,\
                     'gThresholdLocal':1e-8,'gThresholdGlobal':1e-8}

# generate the vortex blobs
x,y = -np.ones(nBlobs)*5.0, np.linspace(-4.0,2.0,nBlobs)
g = np.zeros(nBlobs)

wField = (x,y,g) # the vorticity field made up of the blobs coordinates and circulation

# Generate the blobs
blobs = pHyFlow.vortex.VortexBlobs(wField,vInf,nu,deltaTc,h,overlap,
                                   blobControlParams=blobControlParams)
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Define geometries

# Cylinder 1 (Radius = 1, nPanels = 100)    
R1       = 1.0   # Radius of cylinder
nPanel1  = 100   # Number of panels
dPanel1  = np.spacing(100) # Spacing between panel and colloc. point
theta1   = np.linspace(np.pi,-np.pi,nPanel1+1) # Panel polar angles
dtheta1  = theta1[1]-theta1[0] # Angle spacing
r1       = (R1 + dPanel1) / np.cos(dtheta1/2.0) # Radial location of the panel end points

# Panel Coordinates in cartesian coordinates
xPanel1 = r1*np.cos(theta1 - dtheta1/2)
yPanel1 = r1*np.sin(theta1 - dtheta1/2)

# Make the cylinder and append the parameters to a dictionary.
cylinder1Data = {'xPanel' : r1*np.cos(theta1 - dtheta1/2),
                 'yPanel' : r1*np.sin(theta1 - dtheta1/2),
                 'cmGlobal'   : np.array([0.,1.]),
                 'thetaLocal' : 0.,
                 'dPanel' : np.spacing(100)}

# Cylinder 2 (Radius = 0.5, nPanels = 50)      
R2       = 0.5   # Radius of cylinder
nPanel2  = 50   # Number of panels
dPanel2  = np.spacing(100) # Spacing between panel and colloc. point
theta2   = np.linspace(np.pi,-np.pi,nPanel2+1) # Panel polar angles
dtheta2  = theta2[1]-theta2[0] # Angle spacing
r2       = (R2 + dPanel2) / np.cos(dtheta2/2.0) # Radial location of the panel end points

# Make the cylinder and append the parameters to a dictionary.
cylinder2Data = {'xPanel' : r2*np.cos(theta2 - dtheta2/2),
                 'yPanel' : r2*np.sin(theta2 - dtheta2/2),
                 'cmGlobal'   : np.array([-3.,-1.]),
                 'thetaLocal' : np.pi/2.,
                 'dPanel' : np.spacing(100)}
    
# Cylinder 3 (Radius = 2.0, nPanels = 200)
R3       = 2.0   # Radius of cylinder
nPanel3  = 200   # Number of panels
dPanel3  = np.spacing(100) # Spacing between panel and colloc. point
theta3   = np.linspace(np.pi,-np.pi,nPanel3+1) # Panel polar angles
dtheta3  = theta3[1]-theta3[0] # Angle spacing
r3       = (R3 + dPanel3) / np.cos(dtheta3/2.0) # Radial location of the panel end points

# Make the cylinder and append the parameters to a dictionary.
cylinder3Data = {'xPanel' : r3*np.cos(theta3 - dtheta3/2),
                 'yPanel' : r3*np.sin(theta3 - dtheta3/2),
                 'cmGlobal'   : np.array([3.,-1.]),
                 'thetaLocal' : np.pi,
                 'dPanel' : np.spacing(100)}

# Append all the geometries for a multi-body panel solver

# For now, only a single cylinder
geometries = {'cylinderNormal':cylinder1Data,
              'cylinderSmall' :cylinder2Data,
              'cylinderLarge' :cylinder3Data}

# Initialize the panels
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
                         
x,y  = np.meshgrid(np.linspace(-4,6,30),np.linspace(-4,2,40))
vortexPanel.panels.solve(vInf[0],vInf[1])
vx,vy = vortexPanel.panels.evaluateVelocity(x.flatten(),y.flatten())

#------------------------------------------------------------------------------     

#------------------------------------------------------------------------------
# Plot the convection

py.figure(1)
py.title('Multi-body : Convection of zero-valued blobs')
#py.clf()
py.quiver(x.flat,y.flat,vx+vInf[0],vy+vInf[1],scale=100)
py.scatter(blobs.x,blobs.y,c='g',s=30,edgecolor='none',label='blobs')
[py.plot(x,y,k+'-',label=geoName) for x,y,geoName,k in zip(vortexPanel.panels.xyPanelGlobal[0],
                                                            vortexPanel.panels.xyPanelGlobal[1],
                                                            vortexPanel.panels.geometryKeys,['b','g','k'])]
py.grid(True)
py.axis('scaled')
py.axis([-6,7,-5,4])
    
for i in xrange(2000):  
    if i % 50 == 0:     
        py.scatter(blobs.x,blobs.y,c='g',s=30,edgecolor='none',label='blobs')
        py.draw()
    print i
    vortexPanel.evolve()

#------------------------------------------------------------------------------