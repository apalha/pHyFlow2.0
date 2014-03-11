__doc__ = """

pHyFlow vortexBlobs class population control test.

Description
-----------

A set of uniformly distributed vortex particles are generated with a circulation
equal to sqrt(x*x+y*y). These particle are then subject to population control.
The effects of the two population control parameters are tested.
    

:First added:   2014-01-20  
:Last updated:  2014-01-20     
:Copyright:     Copyright (C) 2014 Artur Palha, **pHyFlow**
:License:       GNU GPL version 3 or any later version
      
"""

"""
 
Reviews:  Added flatten to the definition of the x and y coordinates (apalha, 2014-01-21)
          of blobs                                                  


"""

import pHyFlow
import numpy
import pylab

# define the number of blobs
nBlobs = 10000

# generate the input fields
vInf = numpy.array([1.0, 0.0])              # the free stream velocity
overlap = 1.0                               # the overlap ration between the blobs
h = 1.0/numpy.sqrt(nBlobs)                  # the cell size to which each blob is associated to
                                            # we set h = sqrt(A/nBlobs), where A is the area inside
                                            # which blobs are randomly generated
deltaTc = 0.03                              # the size of the time step, irrelevant for this example as no time stepping is done
nu = 0.001                                  # the dinamic viscous constant, irrelevant for this example as no time stepping is done

# the parameters
# the minimum circulation per blob is set to a small value
R_remove = 0.5
gMin = R_remove*h*h # all blobs located at a place with vorticity smaller than R_remove will be removed
               # given the spatial vorticity defined here, means that a circle of radius R_remove
               # around the origin will be void of particles
gTotalMinVar = 10.0 # we set the total circulation variation to a large value such
                  # that it will not be limiting blob removal

# the parameters
blobControlParams = {'stepRedistribution':1,'stepPopulationControl':1,\
                       'gThresholdLocal':0.5e-2,'gThresholdGlobal':1e-2}

# generate the vortex blobs uniformly distributed in the square [0,1]x[0,1] and
# circulation g = sqrt(x*x+y*y)
x,y = numpy.meshgrid(numpy.linspace(0.0,1.0,numpy.sqrt(nBlobs)),numpy.linspace(0.0,1.0,numpy.sqrt(nBlobs)))
x = x.flatten()
y = y.flatten()
g = h*h*numpy.sqrt(x*x + y*y)
                  
wField = (x,y,g)                            # the vorticity field made up of the blobs coordinates and circulation

# generate the blobs
blobs = pHyFlow.blobs.Blobs(wField,vInf,nu,deltaTc,h,overlap,
                                   blobControlParams=blobControlParams)

# plot the original blobs
pylab.scatter(blobs.x,blobs.y,c=blobs.g,edgecolor='none')

# peform population control
blobs.populationControl()

# plot the redistributed blobs
pylab.figure()
pylab.scatter(blobs.x,blobs.y,c=blobs.g,edgecolor='none')
