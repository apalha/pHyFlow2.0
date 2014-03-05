"""
Test Multi-Body Problem. 3 Cylinders of different size, with different
number of panels.

:First Added:   2013-11-19
"""

import pHyFlow

import numpy as np
import time
import os

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

#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Move the bodies and save the panels in VTK, solving for panel strength every time

# define the filename
os.mkdir('data')
filename = './data/panels.pvd'

# open the pvd file
pvdFile = pHyFlow.IO.File(filename)

# define the time step for the motion functions
dt = 0.1

# step in time
for tStep in range(0,100):
    print '-----------------------------------'
    print 'Time step %d\n' % tStep
    
    startTime_MoveBodies = time.time()
    
    # move the bodies
    cmGlobalNew = {'cylinderNormal': np.array([0.75*np.cos(2*tStep*dt),0.25*np.sin(2*tStep*dt)]),
                   'cylinderSmall' : np.array([3.*np.cos(tStep*dt),3.*np.sin(tStep*dt)]),
                   'cylinderLarge' : np.array([-4.0*np.cos(-tStep*dt),4.0*np.sin(-tStep*dt)])}
    
    # do not rotate them (keep the rotation angle theta constant)
    thetaLocalNew = {'cylinderNormal': 0.,
                     'cylinderSmall' : 0.,
                     'cylinderLarge' : 0.}

    # update panel body
    panelBodies.updateBody(cmGlobalNew,thetaLocalNew)
    panelBodies._advanceTime(dt)
    
    endTime_MoveBodies = time.time()
    
    startTime_UpdateStrengths = time.time()
    
    # update the strengths
    
    # Collocation points
    xCP,yCP = panelBodies.xyCPGlobalCat

    # External velocity at collocation points
    vx,vy = externVel(xCP,yCP)

    # External Velocity : simple free-stream flow
    panelBodies.solve(vx,vy)
    
    endTime_UpdateStrengths = time.time()    
    
    startTime_SaveBodies = time.time()
    
    # save the bodies to pvd
    pvdFile << panelBodies
    
    endTime_SaveBodies = time.time()
    
    # print the times
    print 'Time to move bodies     : %fs' % (endTime_MoveBodies-startTime_MoveBodies)
    print 'Time to update strengths: %fs' % (endTime_UpdateStrengths-startTime_UpdateStrengths)
    print 'Time to save bodies     : %fs' % (endTime_SaveBodies-startTime_SaveBodies)
    print '-----------------------------------\n'