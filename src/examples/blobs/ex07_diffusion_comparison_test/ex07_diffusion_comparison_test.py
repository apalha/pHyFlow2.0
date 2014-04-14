__doc__ = """

pHyFlow vortexBlobs class diffusion test comparison

Description
-----------

A set of uniformly distributed vortex particles are generated with a vorticity
equal to exp(-(x*x + y*y)/(R*R)). This initial vorticity field is then evolved
for a certain amount of time steps and compared to the initial distribution.
Two different diffusion methods are used regrid_wee and regrid_tutty and compared.
    

:First added:   2014-04-14  
:Last updated:  2014-04-14    
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

pylab.ion()

# plot flag that determines if plots are made or not
plot_flag = False
redistribute_flag = True
popcontrol_flag = True
plotanimation_flag = True
plotnumBlobs_flag = True
plotcirculation_flag = True

# define the number of blobs
nBlobs = 64*64#10000
nPlotPoints = 128*128

# generate the input fields
vInf = numpy.array([0.0, 0.0])              # the free stream velocity
overlap = 1.0                               # the overlap ration between the blobs
h = 2.0/numpy.sqrt(nBlobs)                  # the cell size to which each blob is associated to
                                            # we set h = sqrt(A/nBlobs), where A is the area inside
                                            # which blobs are randomly generated
deltaTc = 0.1                               # the size of the time step, irrelevant for this example as no time stepping is done
nu = 0.01                                    # the dinamic viscous constant, irrelevant for this example as no time stepping is done

nTimeSteps = 150                             # evolve the blobs for nTimeSteps


# the parameters
blobControlParams = {'stepRedistribution':1,'stepPopulationControl':1,\
                       'gThresholdLocal':1e-8,'gThresholdGlobal':1e-8}

blobDiffusionParams_tutty = {'method':'regrid_tutty'}
blobDiffusionParams_wee = {'method':'regrid_wee','c2':'optimized'}

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
blobs_tutty = pHyFlow.blobs.Blobs(wField,vInf,nu,deltaTc,h,overlap,
                                  blobControlParams=blobControlParams,diffusionParams=blobDiffusionParams_tutty)

blobs_wee = pHyFlow.blobs.Blobs(wField,vInf,nu,deltaTc,h,overlap,
                                  blobControlParams=blobControlParams,diffusionParams=blobDiffusionParams_wee)
                                  
# change the free stream velocity
blobs_tutty.vInf = numpy.array([0.0,0.0])
blobs_wee.vInf = numpy.array([0.0,0.0])


blobs_pvd_tutty = pHyFlow.IO.File('./results/blobs_tutty.pvd')
blobs_pvd_wee = pHyFlow.IO.File('./results/blobs_wee.pvd')

# evolve the blobs for nTimeSteps
for timeStep in range(0,nTimeSteps):
    print '\n--------------------------------------------'
    print 'Time step tutty :: %d (t=%.4fs)' % (blobs_tutty.tStep+1,blobs_tutty.t+blobs_tutty.deltaTc) # notice that the time information is only updated in the end of the evolve function
    print 'Time step wee :: %d (t=%.4fs)' % (blobs_wee.tStep+1,blobs_wee.t+blobs_tutty.deltaTc) # notice that the time information is only updated in the end of the evolve function
        
    startTime_tutty_evolve = time.time()
    blobs_tutty.evolve() # evolve the blobs
    endTime_tutty_evolve = time.time()
    
    startTime_wee_evolve = time.time()
    blobs_wee.evolve() # evolve the blobs
    endTime_wee_evolve = time.time()
        
    if redistribute_flag:
        startTime_redistribute = time.time()
        blobs_tutty.redistribute() # redistribute the blobs
        blobs_wee.redistribute() # redistribute the blobs
        endTime_redistribute = time.time()
        
    if popcontrol_flag:
        startTime_populationControl = time.time()
        blobs_tutty.populationControl() # perform population control
        blobs_wee.populationControl() # perform population control
        endTime_populationControl = time.time()
        
    

    # display times in each operation of the time step and the number of blobs
    print 'Time to evolve (tutty) : %fs' % (endTime_tutty_evolve - startTime_tutty_evolve)
    print 'Time to evolve (wee)   : %fs' % (endTime_wee_evolve - startTime_wee_evolve)
    
    if redistribute_flag:
        print 'Time to redistribute: %fs' % (endTime_redistribute - startTime_redistribute)
    else:
        print 'Time to redistribute: no redistribution'
    
    if popcontrol_flag:    
        print 'Time to pop control : %fs' % (endTime_populationControl - startTime_populationControl)
    else:
        print 'Time to pop control : no population control'
        
    # display blobs information    
        
    print '\nNumber of blobs     : %d   %d'  % (blobs_tutty.numBlobs, blobs_wee.numBlobs)
    print 'Total circulation     : %.8f   %.8f'  % (blobs_tutty.g.sum(), blobs_wee.g.sum())
    print '--------------------------------------------\n'
        
    # save the blobs to files
    blobs_pvd_tutty << blobs_tutty
    blobs_pvd_wee << blobs_wee

# wait for user to press enter to exit the program
raw_input("Press Enter to continue...")