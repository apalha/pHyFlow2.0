__doc__ = """

pHyFlow vortexBlobs class save time step data to VTK file to open with
Paraview, MayaVI or VTK.

Description
-----------

This example is almost identical to ex04_evolve_test, the only difference is
that blobs are saved to VTK files for plotting.

A set of uniformly distributed vortex particles are generated with a vorticity
equal to exp(-(x*x + y*y)/(R*R)). This initial vorticity field is then evolved
for a certain amount of time steps and compared to the initial distribution.
    

:First added:   2014-02-19  
:Last updated:  2014-02-19    
:Copyright:     Copyright (C) 2014 Artur Palha, **pHyFlow**
:Licence:       GNU GPL version 3 or any later version
      
"""

"""
Reviews:
-------
          

"""

import lxml.etree as lET
from evtk.hl import pointsToVTK
import time
import numpy
import os
import pHyFlow
import pylab

pylab.ion()



# plot flag that determines if plots are made or not
plot_flag = False
redistribute_flag = True
popcontrol_flag = True
plotanimation_flag = True
plotnumBlobs_flag = True
plotcirculation_flag = True
plotVelocityFlag = True

# define the number of blobs
nBlobs = 128*128
nPlotPoints = 128*128

# generate the input fields
vInf = numpy.array([0.0, 0.0])              # the free stream velocity
overlap = 1.0                               # the overlap ration between the blobs
h = 2.0/numpy.sqrt(nBlobs)                  # the cell size to which each blob is associated to
                                            # we set h = sqrt(A/nBlobs), where A is the area inside
                                            # which blobs are randomly generated
deltaTc = 0.01                              # the size of the time step, irrelevant for this example as no time stepping is done
nu = 0.0                                    # the dinamic viscous constant, irrelevant for this example as no time stepping is done

nTimeSteps = 20                             # evolve the blobs for nTimeSteps


# the parameters
blobControlParams = {'stepRedistribution':1,'stepPopulationControl':1,\
                       'gThresholdLocal':1e-8,'gThresholdGlobal':1e-8}
                       
if not redistribute_flag:
    blobControlParams['stepRedistribution'] = 0                

if not popcontrol_flag:
    blobControlParams['stepPopulationControl'] = 0

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
blobs = pHyFlow.vortex.VortexBlobs(wField,vInf,nu,deltaTc,h,overlap,
                                   blobControlParams=blobControlParams)


# change the free stream velocity
blobs.vInf = numpy.array([10.0,0.0])

# change the flag to plot velocity
blobs.plotVelocity = plotVelocityFlag


# if the time evolution of the number of blobs is to be plotted
# allocate memory space for the number of blobs in time
if plotnumBlobs_flag:
    numBlobsTime = numpy.zeros(nTimeSteps+1) # +1 to count for the initial time step
    numBlobsTime[0] = blobs.numBlobs

# if the time evolution of the total circulation is to be plotted
# allocate memory space for the total circulation in time
if plotcirculation_flag:
    totalCirculationTime = numpy.zeros(nTimeSteps+1) # +1 to count for the initial time step
    totalCirculationTime[0] = blobs.g.sum()

filename = './blobs.pvd'
myPVDFile = pHyFlow.IO.File(filename)

myPVDFile << blobs

# evolve the blobs for nTimeSteps
for timeStep in range(0,nTimeSteps):
    print '\n--------------------------------------------'
    print 'Time step :: %d (t=%.4fs)' % (blobs.tStep+1,blobs.t+blobs.deltaTc) # notice that the time information is only updated in the end of the evolve function

    startTime_evolve = time.time()
    blobs.evolve() # evolve the blobs
    endTime_evolve = time.time()

    # get the number of blobs
    if plotnumBlobs_flag:
        numBlobsTime[timeStep+1] = blobs.numBlobs

    # get the total circulation
    if plotcirculation_flag:
        totalCirculationTime[timeStep+1] = blobs.g.sum()
    
    # save the blobs into VTK format    
    startTime_PVDsave = time.time()
    myPVDFile << blobs
    endTime_PVDsave = time.time()
    
    # display times in each operation of the time step and the number of blobs
    print 'Time to evolve      : %fs' % (endTime_evolve - startTime_evolve)

    print 'Time to save to VTK : %fs' % (endTime_PVDsave - startTime_PVDsave)    
    
    # display blobs information

    print '\nNumber of blobs     : %d'  % blobs.numBlobs
    print 'Total circulation     : %.8f'  % blobs.g.sum()
    print '--------------------------------------------\n'

if plotnumBlobs_flag:
    pylab.figure()
    pylab.plot(numpy.arange(0,nTimeSteps+1)*blobs.deltaTc,numBlobsTime)
    pylab.draw()

if plotcirculation_flag:
    pylab.figure()
    pylab.plot(numpy.arange(0,nTimeSteps+1)*blobs.deltaTc,totalCirculationTime-totalCirculationTime[0])
    pylab.draw()
