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
# Define external velocity field

# Free-stream flow
def externVel(x,y):
    vInf = np.array([1.,0.])
    return vInf[0]*np.ones(x.shape[0]), vInf[1]*np.ones(y.shape[0])

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

#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Append all the geometries for a multi-body panel solver

# For now, only a single cylinder
geometries = {'cylinderNormal':cylinder1Data,
              'cylinderSmall' :cylinder2Data,
              'cylinderLarge' :cylinder3Data}

#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Initialize the panels

# Input method 1: Provide all the parameters
# Initalize panelBody
panelBodies = pHyFlow.panels.Panels(geometries=geometries)

#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
# Solve for panel strengths (to satisfy no-slip b.c.)

# Collocation points
xCP,yCP = panelBodies.xyCPGlobalCat

# External velocity at collocation points
vx,vy = externVel(xCP,yCP)

# External Velocity : simple free-stream flow
panelBodies.solve(vx,vy)

py.figure(1)
py.imshow(panelBodies.A,interpolation='none')
py.title('Inter-Induction Matrix')
py.colorbar()
py.show(block=True)

#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Plot all parameters

# Plot the body                  
py.figure(2)
[py.plot(x,y,k+'.-',label=geoName) for x,y,geoName,k in zip(panelBodies.xyPanelGlobal[0],panelBodies.xyPanelGlobal[1],panelBodies.geometryKeys,['b','g','k'])]
py.axis([-10,10,-10,10])
py.grid()
py.legend()
py.axis('scaled')
py.title('Multi-body')
#------------------------------------------------------------------------------
    
#------------------------------------------------------------------------------
# Plot the velocity fields
# Generate grid points
xGrid,yGrid = np.meshgrid(np.linspace(-5,5,100),np.linspace(-3,3,50))

# Calculate induced velocities
vxPanel,vyPanel = panelBodies.evaluateVelocity(xGrid.flatten(), yGrid.flatten())

# Free-stream velocity
vxInf, vyInf = externVel(xGrid.flatten(),yGrid.flatten())

# Total velocity field
vx,vy = vxPanel + vxInf, vyPanel+vyInf

# Plot velocity field
py.figure(2)
py.quiver(xGrid.flatten(),yGrid.flatten(),vx,vy,np.sqrt(vx*vx+vy*vy),scale=10.0,alpha=0.5)
py.colorbar()
py.axis([-5,5,-3,3])
py.show(block=True)

py.figure(4)
py.contourf(xGrid,yGrid,np.sqrt(vx*vx+vy*vy).reshape(xGrid.shape))
py.colorbar()
py.axis('scaled')
py.axis([-5,5,-3,3])
py.title('Mult-body velocity field')
py.show(block=True)

#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Update panel locations
# New positions 

cmGlobalNew = {'cylinderNormal': np.array([0.,0.]),
               'cylinderSmall' : np.array([3.,0.]),
               'cylinderLarge' : np.array([-4.,0.])}

thetaLocalNew = {'cylinderNormal': 0.,
                 'cylinderSmall' : 0.,
                 'cylinderLarge' : 0.}
# update panel body

startTime = time.time()
panelBodies.updateBody(cmGlobalNew,thetaLocalNew)
print "Time to move to new postion: %g seconds." % (time.time() - startTime)


# Plot the body                  
py.figure(5)
#[py.plot(x,y,'ko') for x,y in zip(panelBodies.xCPGlobal,panelBodies.yCPGlobal)]
#[py.plot(x,y,'b.-') for x,y in zip(panelBodies.xPanelGlobal,panelBodies.yPanelGlobal)]
[py.plot(x,y,k+'.-',label=geoName) for x,y,geoName,k in zip(panelBodies.xyPanelGlobal[0],panelBodies.xyPanelGlobal[1],panelBodies.geometryKeys,['b','g','k'])]
py.axis('scaled')
py.axis([-5,5,-3,3])
py.grid()
py.legend()
py.title('Multi-body to new position')

# External velocity at collocation points
vx,vy = externVel(xCP,yCP)

# Solve for panel strengths (to satisfy no-slip b.c.)
panelBodies.solve(vx,vy)

# Plot the velocity fields
                     
# Generate grid points
xGrid,yGrid = np.meshgrid(np.linspace(-5,5,100),np.linspace(-3,3,50))

# Calculate induced velocities
vxPanel,vyPanel = panelBodies.evaluateVelocity(xGrid.flatten(), yGrid.flatten())

# Free-stream velocity
vxInf, vyInf = externVel(xGrid.flatten(),yGrid.flatten())

# Total velocity field
vx,vy = vxPanel + vxInf, vyPanel+vyInf

py.figure(5)
py.contourf(xGrid,yGrid,np.sqrt(vx*vx+vy*vy).reshape(xGrid.shape))
py.colorbar()

# Plot velocity field
py.figure(5)
py.quiver(xGrid.flatten(),yGrid.flatten(),vx,vy,scale=50.0)
py.show(block=True)

#------------------------------------------------------------------------------
#
##------------------------------------------------------------------------------
## Save panel data
##panelBodies.save('./panelData')
##------------------------------------------------------------------------------
#
##------------------------------------------------------------------------------
## Initialize panel problem from save file
##pBody = panels(externVel, initFile='./panelData.npz')
##------------------------------------------------------------------------------    