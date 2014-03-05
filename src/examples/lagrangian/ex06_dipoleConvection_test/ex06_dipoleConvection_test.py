"""

Hybrid test : ex02_dipoleConvection_test
=================================

Goal
----
 - test the coupling of navier-stokes and vortex blobs
     without any panel body.

 - test the evolution of lamb-oseen vortex core     

:First added:   2014-03-03
:Last updated:  2014-03-03
:Copyright:     Copyright (C) 2014 Lento Manickathan, **pHyFlow**
:Licence:       GNU GPL version 3 or any later version
      
"""

# External modules
import os
import numpy as np
import pylab as py
import time
py.ion()


# pHyFlow
import pHyFlow



#-----------------------------------------------------------------------------    
# Global parameters

# Simulation parameters
nTimeSteps = 100
# Fluid Parameters

# Double-monopole parameters
Re = 625. # Reynolds number
vInf = np.array([0.,0.]) # Free stream velocity
U = 1.0
L = 1.0
Tend = 1.0

# Compute fluid kinematic viscosity
nu = U*L/Re         

# Initial conditions for a double-Monopole
# Parameters
x1,y1   = 0.1, 0. # Initial core location. `Positive Core`
x2,y2   = -0.1,0. # Initial core location. `Negative Core`
r0      = 0.1      # Half core radius
omega_e = 299.528385375226# Maximum single core voriticty

# Exact Function
def wExactFunction(x,y):
    rr1 = (x-x1)**2 + (y-y1)**2 # Radius 1
    rr2 = (x-x2)**2 + (y-y2)**2 # Radius 2
    return omega_e*(1- rr1/(r0**2))*np.exp(-rr1/(r0**2)) \
          -omega_e*(1- rr2/(r0**2))*np.exp(-rr2/(r0**2))
                                  
#-----------------------------------------------------------------------------    

#-----------------------------------------------------------------------------
# Setup the vortex blobs - Lagrangian domain

# Define blob parameters
overlap = 1.0               # overlap ratio
nBlobs = 128*128*(2*2)# number of blobs
h = 2.0/np.sqrt(nBlobs)    # blob spacing

# Estimate of delta T
#   with optimized c2 of 1/3.
#   delta T of convection should be a multiple of delta diffusion
nConvectionPerDiffusion = 10.
deltaTc = float(np.round((1.0/nConvectionPerDiffusion) * ( (1./3.)*h*h/nu),decimals=3))

# blob control parameters
blobControlParams = {'stepRedistribution':1,'stepPopulationControl':1,\
                     'gThresholdLocal':1e-8,'gThresholdGlobal':1e-8}


# Generate the vortex blobs uniformly distributed in the square [-10,10] x [-10,10]
x,y = np.meshgrid(np.linspace(-1.0,1.0,np.sqrt(nBlobs)+1),
                  np.linspace(y1-1.0,y1+1.0,np.sqrt(nBlobs)+1))
x = x.flatten()
y = y.flatten()
g = wExactFunction(x,y) * (h*h)

if np.sum(x==0) == 0 or np.sum(y==0) == 0:
    raise ValueError('vortex mesh should coincide with (0.,0.) for remeshing grid')

wField = (x,y,g)

# generate the blobs
blobs = pHyFlow.blobs.Blobs(wField,vInf,nu,deltaTc,h,overlap,
                                   blobControlParams=blobControlParams)

blobs.plotVelocity = True
# This problem does not have any panels
panels = None

# Setup the full vortex problem : blobs and others.
lagrangian = pHyFlow.lagrangian.LagrangianSolver(blobs,panels=None)

#------------------------------------------------------------------------------




#------------------------------------------------------------------------------
# Evolve

# Export the blobs

os.mkdir('data')

blobsFileName = './data/blobs.pvd'
blobsFile = pHyFlow.IO.File(blobsFileName)

# Export initial data
blobsFile << blobs

#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Evolve

# Post-processing parameters
numBlobsTime = np.zeros(nTimeSteps+1)
totalCirculationTime = np.zeros(nTimeSteps+1)
numBlobsTime[0] =  lagrangian.blobs.numBlobs
totalCirculationTime[0] =  lagrangian.blobs.g.sum()


for timeStep in range(1,nTimeSteps+1):
    
    
    
    startTime = time.time()
    # Evolve the coupled navier-stokes and vortex method
    lagrangian.evolve()

    print "\nThe current time\t\t: %g" % (timeStep*lagrangian.deltaT)
    print "Time-step\t\t\t: %g" % timeStep
    print 'Time to evolve\t\t\t: %g'  % (time.time() - startTime)
    
    if timeStep % 10 == 0:
        # Export initial data
        blobsFile << lagrangian.blobs
    
    numBlobsTime[timeStep] =  lagrangian.blobs.numBlobs
    totalCirculationTime[timeStep] =  lagrangian.blobs.g.sum()
    print 'Number of blobs\t\t\t: %d'  % numBlobsTime[timeStep]
    print 'Total circulation\t\t: %.8f'  % totalCirculationTime[timeStep]
    print '--------------------------------------------\n'

#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
# Plot results

py.figure()
py.plot(np.arange(0,nTimeSteps+1)*blobs.deltaTc,numBlobsTime)
py.grid()
py.xlabel('t')

py.figure()
py.semilogy(np.arange(0,nTimeSteps+1)*lagrangian.deltaT,np.abs(totalCirculationTime-totalCirculationTime[0]))
py.grid()
py.xlabel('t')

