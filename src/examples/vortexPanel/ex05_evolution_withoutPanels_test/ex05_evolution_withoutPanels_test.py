__doc__ = """

pHyFlow vortexBlobs class evolve blobs test.

Description
-----------

A set of uniformly distributed vortex particles are generated with a vorticity
equal to exp(-(x*x + y*y)/(R*R)). This initial vorticity field is then evolved
for a certain amount of time steps and compared to the initial distribution.
    

:First added:   2014-02-10  
:Last updated:  2014-02-10    
:Copyright:     Copyright (C) 2014 Artur Palha, **pHyFlow**
:Licence:       GNU GPL version 3 or any later version
      
"""

"""
Reviews:
-------
          

"""

import pHyFlow
import numpy as np
import time
import pylab as py

py.ion()

# plot flag that determines if plots are made or not
plot_flag = False
redistribute_flag = True
popcontrol_flag = True
plotanimation_flag = True
plotnumBlobs_flag = True
plotcirculation_flag = True

#------------------------------------------------------------------------------
# Initialize vortex blobs

# define the number of blobs
nBlobs = 128*128#10000
nPlotPoints = 128*128

# generate the input fields
vInf = np.array([0.0, 0.0])              # the free stream velocity
overlap = 1.0                               # the overlap ration between the blobs
h = 2.0/np.sqrt(nBlobs)                  # the cell size to which each blob is associated to
                                            # we set h = sqrt(A/nBlobs), where A is the area inside
                                            # which blobs are randomly generated
deltaTc = 0.001                              # the size of the time step, irrelevant for this example as no time stepping is done
nu = 0.01#0.005                                  # the dinamic viscous constant, irrelevant for this example as no time stepping is done

nTimeSteps = 100                             # evolve the blobs for nTimeSteps


# the parameters
blobControlParams = {'stepRedistribution':1,'stepPopulationControl':1,\
                     'gThresholdLocal':1e-8,'gThresholdGlobal':1e-8}

# default parameters are used

print 'Generating blob definitions...'

# generate the vortex blobs uniformly distributed in the square [-1,1]x[-1,1] and
# and with associated vorticity omega = exp(-(x*x + y*y)/(R*R))
R = 0.2 # the radius of the Gaussian distribution of vorticity
x,y = np.meshgrid(np.linspace(-1.0,1.0,np.sqrt(nBlobs)),np.linspace(-1.0,1.0,np.sqrt(nBlobs)))
x = x.flatten()
y = y.flatten()
g = h*h*100.0*np.exp(-(x*x + y*y)/(R*R))
                  
wField = (x,y,g)                            # the vorticity field made up of the blobs coordinates and circulation

print 'Generating blobs object...'

# generate the blobs
blobs = pHyFlow.vortex.VortexBlobs(wField,vInf,nu,deltaTc,h,overlap,
                                   blobControlParams=blobControlParams)


# change the free stream velocity
blobs.vInf = np.array([50.0,0.0])

#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Initialize the vortex-panel class
                       
startTime = time.time()                     
vortexPanel = pHyFlow.vortexPanel.VortexPanel(blobs,panels=None,
                                              couplingParams={'panelStrengthUpdate':'varying'})
print "\nTime to initialize the coupled vortex-panel problem: %g seconds." % (time.time() - startTime)
                         
                         
#------------------------------------------------------------------------------
                         
#------------------------------------------------------------------------------
# Plot the convection

py.figure(1)

    
for i in xrange(nTimeSteps): 
    
    if i % 5 == 0:
        py.clf()
        py.title('Convection of gaussian field')
        py.scatter(vortexPanel.blobs.x,vortexPanel.blobs.y,c=vortexPanel.blobs.g,edgecolor='none')
        py.grid(True)
        py.axis('scaled')
        py.axis([-1,5,-1.5,1.5])
        py.draw()
    print i
    
    # Convect only the blobs
    vortexPanel.evolve()

#------------------------------------------------------------------------------             