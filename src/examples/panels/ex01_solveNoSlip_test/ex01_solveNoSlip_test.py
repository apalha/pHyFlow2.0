
__doc__=r"""
Testing single cylinder case. Compare the numerical result with
analytical results.
"""

import pHyFlow
import numpy as np
import pylab as py
py.ion()
import time



#------------------------------------------------------------------------------
# Free-stream flow
def externVel(x,y):
    vInf = np.array([1.,0.])
    return vInf[0]*np.ones(x.shape[0]), vInf[1]*np.ones(y.shape[0])  
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
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Append all the geometries for a multi-body panel solver

# For now, only a single cylinder
geometries = {'cylinder':cylinderData}

#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Initialize the panels

startTime = time.time()

# Initalize panelBody 
panelBodies = pHyFlow.panels.Panels(geometries=geometries)

print "\nTime to initialize the panel problem: %g seconds." % (time.time() - startTime)

# Self induction Matrix
A = panelBodies.A

py.figure(1)
py.imshow(A,interpolation='none')
py.title('Self-Induction Matrix (%g by %g)' % (nPanel, nPanel))
py.colorbar()
py.show(block=True)

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
# Solve for the panel strengths

# Panel Strength
sPanel = panelBodies.sPanel[0]

py.figure(2)
py.plot(sPanel,'k.-')
py.title('Panel Strengths')
py.xlabel('Panel Number')
#------------------------------------------------------------------------------
