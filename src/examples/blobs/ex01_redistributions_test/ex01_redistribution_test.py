__doc__ = """

pHyFlow vortexBlobs class redistribution test.

Description
-----------

A set of randomly distributed vortex particles are generated with a circulation
equal to sin(2*pi*x)*sin(2*pi*y). These particle are then redistributed into a
regular grid using the redistribute method of vortexBlobs class.
    

:First added:   2014-01-16  
:Last updated:  2014-01-16     
:Copyright:     Copyright (C) 2014 Artur Palha, **pHyFlow**
:Licence:       GNU GPL version 3 or any later version
      
"""

"""
Reviews:
-------
          

"""

import pHyFlow
import numpy
import pylab

# define the number of blobs
nBlobs = 10000
perturbation = (1.0/numpy.sqrt(nBlobs))*1.0

# generate the vortex blobs uniformly distributed in the square [0,1]x[0,1] and
# circulation g = sin(2*pi*x)*sin(2*pi*y)
x,y = numpy.meshgrid(numpy.linspace(0.0,1.0,numpy.sqrt(nBlobs)),numpy.linspace(0.0,1.0,numpy.sqrt(nBlobs)))
x = x.flatten() + perturbation * (numpy.random.rand(nBlobs)-0.5)/0.5
y = y.flatten() + perturbation * (numpy.random.rand(nBlobs)-0.5)/0.5
g = numpy.sin(2.0*numpy.pi*x)*numpy.sin(2.0*numpy.pi*y)*(1.0/numpy.sqrt(nBlobs))*(1.0/numpy.sqrt(nBlobs))

# add a small random perturbation to the positions of all particles

# plot the randomly distributed blobs
pylab.scatter(x,y,c=g,edgecolor='none')

# generate the input fields
wField = (x,y,g)                            # the vorticity field made up of the blobs coordinates and circulation
vInf = numpy.array([1.0, 0.0])              # the free stream velocity
overlap = 1.0                               # the overlap ration between the blobs
h = 1.0/numpy.sqrt(nBlobs)                  # the cell size to which each blob is associated to
                                            # we set h = sqrt(A/nBlobs), where A is the area inside
                                            # which blobs are randomly generated
deltaTc = 0.03                              # the size of the time step, irrelevant for this example as no time stepping is done
nu = 0.001                                  # the dinamic viscous constant, irrelevant for this example as no time stepping is done

# the parameters are not specified therefore default_params are used

# generate the blobs
blobs = pHyFlow.blobs.Blobs(wField,vInf,nu,deltaTc,h,overlap)

# redistribute the blobs
blobs.redistribute()

# plot the redistributed blobs
pylab.figure()
pylab.scatter(blobs.x,blobs.y,c=blobs.g,edgecolor='none')
