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
import numpy
import pylab
import time


# plot flag that determines if plots are made or not
plot_flag = True


# define the number of blobs
nBlobs = 128*128#10000
nPlotPoints = 128*128

# generate the input fields
vInf = numpy.array([0.0, 0.0])              # the free stream velocity
overlap = 1.0                               # the overlap ration between the blobs
h = 2.0/numpy.sqrt(nBlobs)                  # the cell size to which each blob is associated to
                                            # we set h = sqrt(A/nBlobs), where A is the area inside
                                            # which blobs are randomly generated
deltaTc = 0.01                              # the size of the time step, irrelevant for this example as no time stepping is done
nu = 0.0                                    # the dinamic viscous constant, irrelevant for this example as no time stepping is done



# default parameters are used

print 'Generating blob definitions...'

# generate the vortex blobs uniformly distributed in the square [-1,1]x[-1,1] and
# and with associated vorticity omega = exp(-(x*x + y*y)/(R*R))
R = 0.2 # the radius of the Gaussian distribution of vorticity
x,y = numpy.meshgrid(numpy.linspace(-1.0,1.0,numpy.sqrt(nBlobs)),numpy.linspace(-1.0,1.0,numpy.sqrt(nBlobs)))
x = x.flatten()
y = y.flatten()
g = h*h*100.0*numpy.exp(-(x*x + y*y)/(R*R))
                  
wField = (x,y,g)                            # the vorticity field made up of the blobs coordinates and circulation

print 'Generating blobs object...'

# generate the blobs
blobs = pHyFlow.blobs.Blobs(wField,vInf,nu,deltaTc,h,overlap)

# get the original number of blobs
numBlobsOriginal = blobs.numBlobs

# plot the original blobs
if plot_flag:
    pylab.figure()    
    pylab.scatter(blobs.x,blobs.y,c=blobs.g,edgecolor='none')
    pylab.axis('equal')    
    pylab.xlim([-1.0,1.0])
    pylab.ylim([-1.0,1.0])

# remove all the blobs that are on the positive half plain (x>0)

# generate the bool iBlobs selection array
iBlobs = (blobs.x < 0)

# remove the blobs
startTimeRemove = time.time()
blobs.removeBlobs(iBlobs)
endTimeRemove = time.time()

# get the number of blobs after removal
numBlobsAfterRemoval = blobs.numBlobs

# plot the new blobs
if plot_flag:
    pylab.figure()    
    pylab.scatter(blobs.x,blobs.y,c=blobs.g,edgecolor='none')
    pylab.axis('equal')    
    pylab.xlim([-1.0,1.0])
    pylab.ylim([-1.0,1.0])
    
# display information

print '\n---------------------------------------'
print 'Number of blobs original       : %d' % numBlobsOriginal
print 'Number of blobs after removal  : %d' % numBlobsAfterRemoval
print '\nTime of removal              : %.8f' % (endTimeRemove-startTimeRemove)
print '-----------------------------'

