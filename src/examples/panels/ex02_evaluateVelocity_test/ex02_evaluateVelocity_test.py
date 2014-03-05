"""
Testing single cylinder case. Compare the numerical result with
analytical results.
"""


import pHyFlow
import time
import numpy as np
import pylab as py
py.ion()





#------------------------------------------------------------------------------
# Free-stream flow
def externVel(x,y):
    vInf = np.array([1.,0.])
    return vInf[0]*np.ones(x.shape[0]), vInf[1]*np.ones(y.shape[0])  
#------------------------------------------------------------------------------
    
#------------------------------------------------------------------------------
# Define the panel geometries : cylinder
    
# Cylinder (Radius = 1, nPanels = 100)    
R       = 1.0   # Radius of cylinder
nPanel  = 100   # Number of panels
dPanel  = np.spacing(100) # Spacing between panel and colloc. point
theta  = np.linspace(np.pi,-np.pi,nPanel+1) # Panel polar angles
dtheta      = theta[1]-theta[0] # Angle spacing
r           = (R + dPanel) / np.cos(dtheta/2.0) # Radial location of the panel end points

# Make the cylinder and append the parameters to a dictionary.
cylinderData = {'xPanel' : r*np.cos(theta - dtheta/2),
                'yPanel' : r*np.sin(theta - dtheta/2),
                'cmGlobal'   : np.array([0.,0.]),
                'thetaLocal' : 0.,
                'dPanel' : np.spacing(100)}
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Append all the geometries for a multi-body panel solver

# For now, only a single cylinder
geometries = {'cylinder':cylinderData}

#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Initialize the panels

# Initalize panelBody 
panelBodies = pHyFlow.panels.Panels(geometries=geometries)

#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Solve for the panel strengths

# Collocation points
xCP,yCP = panelBodies.xyCPGlobalCat

# External velocity at collocation points
vx,vy = externVel(xCP,yCP)

startTime = time.time()
panelBodies.solve(vx,vy)
print "\nTime to solve for panel strengths: %g seconds." % (time.time() - startTime)

#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Plot the velocity fields

# Generate grid points
xGrid,yGrid = np.meshgrid(np.linspace(-2,2,100),np.linspace(-2,2,50))

# Calculate induced velocities
startTime = time.time()
vxPanel,vyPanel = panelBodies.evaluateVelocity(xGrid.flatten(), yGrid.flatten())
print "\nTime to evaluate the velocity field: %g seconds." % (time.time() - startTime)
# Free-stream velocity
vxInf, vyInf = externVel(xGrid.flatten(),yGrid.flatten())

# Total velocity field
vx,vy = vxPanel+vxInf, vyPanel+vyInf
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Retrive geometry data

# Panel Global collocation points: in list
xCPGlobal, yCPGlobal = panelBodies.xyCPGlobal

# Panel corner points: in list
xPanelGlobal, yPanelGlobal = panelBodies.xyPanelGlobal

# Geometry keys
geometryKeys = panelBodies.geometryKeys
#------------------------------------------------------------------------------

#-------------------------------1-----------------------------------------------
# Plots

# Plot velocity field
py.figure(1)
[py.plot(x,y,'r.',label=geoName + ' [colloc. points]') for x,y,geoName in zip(xCPGlobal,yCPGlobal,geometryKeys)]
[py.plot(x,y,'k.-',label=geoName + ' [panel]') for x,y,geoName in zip(xPanelGlobal,yPanelGlobal,geometryKeys)]
py.grid()
py.legend()
py.quiver(xGrid.flatten(),yGrid.flatten(),vx,vy,np.sqrt(vx*vx+vy*vy),scale=10.0,alpha=0.5)
py.axis('scaled'), py.title('Cylinder'), py.colorbar(), py.axis([-2,2,-2,2])
py.show(block=True)

py.figure(2)
[py.plot(x,y,'r.',label=geoName + ' [colloc. points]') for x,y,geoName in zip(xCPGlobal,yCPGlobal,geometryKeys)]
[py.plot(x,y,'k.-',label=geoName + ' [panel]') for x,y,geoName in zip(xPanelGlobal,yPanelGlobal,geometryKeys)]
py.contourf(xGrid,yGrid,np.sqrt(vx*vx+vy*vy).reshape(xGrid.shape))
py.axis('scaled')
py.axis([-2,2,-2,2])
py.colorbar()
py.title('Cylinder velocity field')
py.show(block=True)
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Calculate the error : w.r.t to the analytical solution : along the y-axis from cylinder boundary.

# Calculate the induced velocity along the y-axis from cylinder boundary.
nPoints = 300
Rinf    = 10.0
#xEval, yEval = 0.0 * np.ones(nPoints), np.linspace(R+0.01,Rinf,nPoints)
xEval, yEval = 0.0 * np.ones(nPoints), np.linspace(R,Rinf,nPoints)

# Calculate induced velocities on the y-axis
vxPanel,vyPanel = panelBodies.evaluateVelocity(xEval, yEval)

# Free-stream velocity
vxInf, vyInf = externVel(xEval,yEval)

# Total velocity field
vx,vy = vxPanel + vxInf, vyPanel+vyInf

# Analytical equation
vxAnalytical = externVel(xEval,yEval)[0]*( 1.0 + ((R**2)/(yEval**2)))

# Plot the solutions: analytical vs. numerical
py.figure(3)
py.title('Velocity along y-axis')
py.plot(yEval,vx,'b.-',label='numerical')
py.plot(yEval,vxAnalytical,'r.-',label='analytical')
py.grid(), py.legend(loc=0)
py.show(block=True)

# Plot the error: analytical vs. numerical : logscale
py.figure(4)
py.title('Error in Velocity along y-axis')
py.semilogy(yEval, np.abs(vx - vxAnalytical)/np.abs(vxAnalytical),'b.-')
py.grid()
py.show(block=True)
#------------------------------------------------------------------------------    
