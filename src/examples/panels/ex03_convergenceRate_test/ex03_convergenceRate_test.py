"""
Testing single cylinder case. Compare the numerical result with
analytical results.
"""

import pHyFlow
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
    
# For testing convergence    
def get_inducedVelocity_of_cylinder(nPanel,xEval,yEval):
    """
    Cylinder : 
        - Radius = 1
        - cmGlobal = [0.,0.]
        - thetaLocal = 0.
    """
    
    # Define the geometyr
    theta   = np.linspace(np.pi,-np.pi,nPanel+1) # Panel polar angles
    dtheta  = theta[1]-theta[0] # Angle spacing
    r       = (R + dPanel) / np.cos(dtheta/2.0) # Radial location of the panel end points
  
    # Make the cylinder and append the parameters to a dictionary.
    cylinderData = {'xPanel' : r*np.cos(theta - dtheta/2),
                    'yPanel' : r*np.sin(theta - dtheta/2),
                    'cmGlobal'   : np.array([0.,0.]),
                    'thetaLocal' : 0.,
                    'dPanel' : np.spacing(100)}
    # Initialize panelBody 
    panelBodies = pHyFlow.panels.Panels(geometries={'cylinder': cylinderData})
    
    # Free-stream flow
    vxInf,vyInf = externVel(panelBodies.xyCPGlobalCat[0],
                            panelBodies.xyCPGlobalCat[1])
    
    # Solve panel body problem
    panelBodies.solve(vxInf,vyInf)
    
    # Calculate induced velocities
    vxPanel,vyPanel = panelBodies.evaluateVelocity(xEval, yEval)
    
    # Free-stream flow
    vxInf,vyInf = externVel(xEval, yEval)
    
    return vxPanel+vxInf,vyPanel+vyInf
#------------------------------------------------------------------------------



#------------------------------------------------------------------------------
# Determing the convergence rate

R       = 1.0   # Radius of cylinder
dPanel  = np.spacing(100) # Spacing between panel and colloc. point

# Number of panels
nPanels = np.array([10,50,100,200,500,1000,2000])    

# Store Data
vxNumerical = np.zeros(7)

# Error evaluation location
xEval, yEval = np.array([0.]), np.array([1.5])

# Calculate the numerical velocity
for i,nPanel in enumerate(nPanels):
    vxNumerical[i],garbage = get_inducedVelocity_of_cylinder(nPanel,xEval,yEval)
    # Calculate analytical velocity        
vxAnalytical = externVel(xEval,yEval)[0]*( 1.0 + ((R**2)/(yEval**2)))

# Plot Error
py.figure()
py.loglog(nPanels,np.abs(vxNumerical-vxAnalytical)/np.abs(vxAnalytical),'b.-')
py.grid()
py.xlabel(r'$N_{panels}$')
py.ylabel(r'Relative Error $\epsilon_r$')
py.title('Convergence of constant strength straight vortex panels')
py.show(block=True)
    
#------------------------------------------------------------------------------



