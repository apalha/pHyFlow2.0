__doc__ = """

pHyFlow vortexBlobs class evaluate velocity test.

Description
-----------

A set of uniformly distributed vortex particles are generated with a vorticity
equal to exp(-(x*x + y*y)/(R*R)). The velocity field generated by this vorticity
distribution is computed and plotted.
    

:First added:   2014-01-21  
:Last updated:  2014-01-21    
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
plot_flag = False

# define the number of blobs
nBlobs = 256*256#10000
nPlotPoints = 256*256

# generate the input fields
vInf = numpy.array([1.0, 0.0])              # the free stream velocity
overlap = 1.0                               # the overlap ration between the blobs
h = 2.0/numpy.sqrt(nBlobs)                  # the cell size to which each blob is associated to
                                            # we set h = sqrt(A/nBlobs), where A is the area inside
                                            # which blobs are randomly generated
deltaTc = 0.0001                            # the size of the time step, irrelevant for this example as no time stepping is done
nu = 0.001                                  # the dinamic viscous constant, irrelevant for this example as no time stepping is done

# the parameters

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
blobs = pHyFlow.vortex.VortexBlobs(wField,vInf,nu,deltaTc,h,overlap)



# plot the original blobs
if plot_flag:
    pylab.figure()    
    pylab.scatter(blobs.x,blobs.y,c=blobs.g,edgecolor='none')

print 'Applying population control...'

## peform population control
blobs.populationControl()

# plot the blobs after population control
if plot_flag:
    pylab.figure()
    pylab.scatter(blobs.x,blobs.y,c=blobs.g,edgecolor='none')


# compute the points where to evaluate velocity to plot it
xPlot,yPlot = numpy.meshgrid(numpy.linspace(-1.0,1.0,numpy.sqrt(nPlotPoints)),numpy.linspace(-1.0,1.0,numpy.sqrt(nPlotPoints)))
xPlot = xPlot.flatten()
yPlot = yPlot.flatten()

# CPU calculation

print 'Computing induced velocity with CPU...'

# set the hardware to cpu
blobs.hardware = 'cpu'

# the first evaluation is done to make more correct evaluation of the time
# since first time the library is ran, the library is loaded and this takes
# time, which is an overhead that happens only once
vxCPU,vyCPU = blobs.evaluateVelocity(xPlot.copy(),yPlot.copy())

# time it
cpu_start = time.time()

# compute the velocity field at the blob points, we can use x,y or blobs.x, blobs.y
vxCPU,vyCPU = blobs.evaluateVelocity(xPlot.copy(),yPlot.copy())

# time it
cpu_time = time.time() - cpu_start

# GPU calculation

print 'Computing induced velocity with GPU...'

# set the hardware to gpu
blobs.hardware = 'gpu'

# the first evaluation is done to make more correct evaluation of the time
# since first time the library is ran, the library is loaded and this takes
# time, which is an overhead that happens only once
vxGPU,vyGPU = blobs.evaluateVelocity(xPlot.copy(),yPlot.copy())

# time it
gpu_start = time.time()

# compute the velocity field at the blob points, we can use x,y or blobs.x, blobs.y
vxGPU,vyGPU = blobs.evaluateVelocity(xPlot.copy(),yPlot.copy())

# time it
gpu_time = time.time() - gpu_start

# plot the velocity field
if plot_flag:
    # plot GPU velocities    
    pylab.figure()
    pylab.quiver(xPlot,yPlot,vxGPU,vyGPU,numpy.sqrt(vxGPU*vxGPU+vyGPU*vyGPU))
    
    # plot CPU velocities
    pylab.figure()
    pylab.quiver(xPlot,yPlot,vxCPU,vyCPU,numpy.sqrt(vxCPU*vxCPU+vyCPU*vyCPU))
    
    # plot the difference in velocities between CPU and GPU
    pylab.figure()
    pylab.quiver(xPlot,yPlot,vxCPU-vxGPU,vyCPU-vyGPU,numpy.sqrt((vxCPU-vxGPU)*(vxCPU-vxGPU)+(vyCPU-vyGPU)*(vyCPU-vyGPU)))
    pylab.colorbar()

# print maximum difference between CPU and GPU calculations
print 'Max difference in x direction: %.5e.' % numpy.abs(vxCPU-vxGPU).max()
print 'Max difference in y direction: %.5e.' % numpy.abs(vyCPU-vyGPU).max()
    
# display the times
print 'CPU time : %.6fs' % cpu_time
print 'GPU time : %.6fs' % gpu_time # 0.5s for nBlobs = 160 000, MATLAB 0.37s
