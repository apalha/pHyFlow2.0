"""

    ex_06
    
        Example 06 showing the usage of time evolution routines for leap frogging
        
        Generate (nBlobs x nBlobs) vortex blobs evenly spaced in the interval 
        [-1, 1] x [-1, 1] with associated circulation given by:
        
            w = gamma*(exp(-(((x-pVortex(1)).^2) + ((y-pVortex(2)).^2))/(2*epsilon*epsilon)))/(2*pi*epsilon*epsilon)
            
            w(x,y) = exp(-(((x+0.1)^2)+y^2)/(C0))/(pi*C0) -
                     -exp(-(((x+0.02)^2)+y^2)/(C0))/(pi*C0) +
                     +exp(-(((x-0.02)^2)+y^2)/(C0))/(pi*C0) -
                     -exp(-(((x-0.1)^2)+y^2)/(C0))/(pi*C0)
                                 
        with:
            C0 = 4*mu*t/rho
            mu = 1.75E-5
            rho = 1.2
            t = 1            
            deltaX = deltaY = 2.0/nBlobs
        
        Advance in time the vortex blobs without and with remeshing.

"""
# Copyright (C) 2013 Artur Palha                                                                                                     
#                                                                                                                                   
# This file is part of pHyFlow.                                                                                                      
#                                                                                                                                   
# pHyFlow is free software: you can redistribute it and/or modify                                                                    
# it under the terms of the GNU Lesser General Public License as published by                                                       
# the Free Software Foundation, either version 3 of the License, or                                                                 
# (at your option) any later version.                                                                                               
#                                                                                                                                   
# pHyFlow is distributed in the hope that it will be useful,                                                                         
# but WITHOUT ANY WARRANTY; without even the implied warranty of                                                                    
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the                                                                      
# GNU Lesser General Public License for more details.                                                                               
#                                                                                                                                   
# You should have received a copy of the GNU Lesser General Public License                                                          
# along with pHyFlow. If not, see <http://www.gnu.org/licenses/>.                                                                    
#                                                                                                                                   
# First added:  2013-07-16                                                                                                          

"""    
    Reviews:
        
"""

import pHyFlow
import numpy
import time
import pylab

#=============================================================================
# input parameters

# parameters definition
nBlobs = numpy.array([256,256])         # the number of vortex blobs in x and y

overlap = 0.5       # the overlap between vortex blobs (together with the number
                    # of particles defines the core size of the particles)

k = 2.0;            # the peakedness of the kernel typical value is 2.0

blockSize = 128      # the size of the block of the GPU, this value is strongly
                    # related to the speed of the calculations on the GPU, for
                    # a large number of particles (> 1024) 128 is a good number.

gthreshold = 10**-5 # the minimum value of circulation to consider
gmaxthreshold = 10**-4 # the maximum value of total circulation variation to consider

dt = 0.01           # the time step size
tSteps = 500        # the number of time steps

stepDisplay = 10    # the interval in time steps between plots of the solution  
stepRedistribute = 10 # the interval in time steps between redistribution

xBoundsDomain = numpy.array([-1.0,1.0]) # define the bounds for the domain
yBoundsDomain = numpy.array([-1.0,1.0]) # define the bounds for the domain

xBoundsPlot = numpy.array([-1.0,3.0]) # define the bounds for the plot
yBoundsPlot = numpy.array([-1.0,1.0]) # define the bounds for the plot

nPlotPoints = [128,128]

mu = 1.75E-5 # viscosity for computing vortex spread
rho = 1.2    # density for computing vortex spread
t = 1.0      # time for computing vortex spread
C0 = 4*mu*t/rho # C0, core radius vortex

def wExact(x,y): return (1.0 * numpy.exp(-(((x+0.1)**2)+(y)**2)/(C0))/(numpy.pi*C0)\
                         -1.0 * numpy.exp(-(((x+0.02)**2)+(y)**2)/(C0))/(numpy.pi*C0)\
                         +1.0 * numpy.exp(-(((x-0.02)**2)+(y)**2)/(C0))/(numpy.pi*C0)\
                         -1.0 * numpy.exp(-(((x-0.1)**2)+(y)**2)/(C0))/(numpy.pi*C0))
#=============================================================================



#=============================================================================
# computation of parameters

h = (xBoundsDomain[1]-xBoundsDomain[0])/(nBlobs[0])   # spacing in x of vortex blobs (note that particles
                                                       # have the same spacing in y => deltaY = deltaX)

sigma = (h/overlap) # compute the core spreading
sigmasqr = sigma * sigma # compute the square of the core spreading
#=============================================================================



#=============================================================================
# computation of vortex blobs

# generate the coordinates of the vortex blobs
xB = xBoundsDomain[0]+0.5*h + numpy.arange(0,nBlobs[0])*h        # nBlobs evenly spaced in x \in [xBoundsDomain[0], xBoundsDomain[1]]
yB = yBoundsDomain[0]+0.5*h + numpy.arange(0,nBlobs[1])*h        # nBlobs evenly spaced in y \in [yBoundsDomain[0], yBoundsDomain[1]]
xB,yB = numpy.meshgrid(xB,yB)                               # generate the 2d grid of blobs in [xBoundsDomain[0], xBoundsDomain[1]] x [yBoundsDomain[0], yBoundsDomain[1]]

# flatten the grids into vectors
xB = xB.flatten()                           
yB = yB.flatten()

# compute the circulations of the blobs
wB = wExact(xB,yB)*h*h
#=============================================================================

# plot vorticity from circulation using dolfin
#figureHandleCirculationAll = pHyFlow.vortex.plotvorticity(xBoundsPlot,yBoundsPlot,xB,yB,wB,\
#                                            sigma,overlap,k=k,\
#                                            wType=pHyFlow.options.WBLOB_VORTICITY,\
#                                            figureHandle=None,title='Vorticity from circulation no pop control',plotType=pHyFlow.options.DOLFIN_PLOT,interactive=True)
                                            
xB,yB,wB = pHyFlow.vortex.regrid.PopulationControl(xB,yB,wB,1e-8,1e-6)

#figureHandleCirculation = pHyFlow.vortex.plotvorticity(xBoundsPlot,yBoundsPlot,xB,yB,wB,\
#                                            sigma,overlap,k=k,\
#                                            wType=pHyFlow.options.WBLOB_VORTICITY,\
#                                            figureHandle=None,title='Vorticity from circulation pop control',plotType=pHyFlow.options.DOLFIN_PLOT,interactive=True)                                                   

# plot the blob centers
pylab.scatter(xB,yB,c=wB)
pylab.axis('equal')    
pylab.axis(numpy.concatenate((xBoundsPlot,yBoundsPlot)))
pylab.draw()

# allocate memory space for the maximum value of vorticity and total circulation
wMax = numpy.zeros(tSteps+1)
wTotal = numpy.zeros(tSteps+1)
nBlobsTime = numpy.zeros(tSteps+1)

# compute the initial maximum value of vorticity and total circulation
wMax[0] = (numpy.abs(wB)).max()/(h*h)
wTotal[0] = wB.sum()

# get the initial number of particles
nBlobsTime[0] = xB.size

for tStep in range(1,tSteps+1):
    print '----------------------------------------------------------'
    print '   Time step :: ' + str(tStep)
    
    start = time.time()
    xB,yB,wB = pHyFlow.vortex.evolution.convection(dt,xB,yB,wB,sigma,k=k,kernel=pHyFlow.options.GAUSS_KERNEL,vInf=[1.0,0.0]\
                                                   hardware=pHyFlow.options.GPU_HARDWARE,method=pHyFlow.options.DIRECT_METHOD,\
                                                   blocksize=128,integrator=pHyFlow.options.RK4_INTEGRATOR)
    print '   Time to convect            :: ' + '%14.10f' % (time.time()-start) + 's'
    
    if (tStep % stepRedistribute) == 0:
        start = time.time()                                         
        xB,yB,wB = pHyFlow.vortex.regrid.Regrid(xB,yB,wB,sigma,overlap,xBoundsDomain,yBoundsDomain,interpKernel=0)
        print '   Time to redistribute       :: ' + '%14.10f' % (time.time()-start) + 's'
    
        xB,yB,wB = pHyFlow.vortex.regrid.PopulationControl(xB,yB,wB,gthreshold,gmaxthreshold)    
    
    
    # plot vorticity from vorticity function using dolfin
    if (tStep % stepDisplay) == 0:    
        #figureHandleCirculation = pHyFlow.vortex.plotvorticity(xBoundsPlot,yBoundsPlot,xB,yB,wB,\
        #                                                 sigma,overlap,k=k,kernel=pHyFlow.options.GAUSS_KERNEL,\
        #                                                 nPlotPoints=nPlotPoints,hardware=pHyFlow.options.GPU_HARDWARE,\
        #                                                 blocksize=128,method=pHyFlow.options.DIRECT_METHOD,\
        #                                                 wType=pHyFlow.options.WBLOB_VORTICITY,\
        #                                                 figureHandle=figureHandleCirculation,title='Vorticity from blob function',plotType=pHyFlow.options.DOLFIN_PLOT,interactive=True) 
        
        pHyFlow.vortex.plotblobs(xBoundsPlot,yBoundsPlot,xB,yB,wB,figureHandle=1,title='Vortex blobs (t=%10.4fs)' % (dt*tStep))
        
    
    wMax[tStep] = ((numpy.abs(wB)).max()/(h*h))
    print '   Maximum value of vorticity :: ' + '%14.10f (%14.10f)' % (wMax[tStep],wMax[tStep]-wMax[0])
    
    wTotal[tStep] = wB.sum()
    print '   Total circulaiton          :: ' + '%14.10f (%14.10f)' % (wTotal[tStep],wTotal[tStep]-wTotal[0])
    
    nBlobsTime[tStep] = xB.size
    print '   Number of particles        :: ' + '%6.0f' % nBlobsTime[tStep]
    
    
    print '----------------------------------------------------------'
    
    
    
# plot the summary of the results
pylab.figure()
pylab.plot(dt*numpy.arange(1,tSteps+1),wMax[1:] - wMax[0],'-ob')
pylab.xlabel(r't [s]')
pylab.ylabel(r'$|\omega_{max}(t)-\omega_{max}(0)|$')
pylab.title('Variation in maximum vorticity')

pylab.figure()
pylab.semilogy(dt*numpy.arange(1,tSteps+1),numpy.abs(wTotal[1:] - wTotal[0]),'-ob')
pylab.xlabel(r't [s]')
pylab.ylabel(r'$|\Gamma_{total}(t)-\Gamma_{total}(0)|$')
pylab.title('Variation in total circulation')

pylab.figure()
pylab.plot(dt*numpy.arange(0,tSteps+1),nBlobsTime,'-ob')
pylab.xlabel(r't [s]')
pylab.ylabel(r'$\sharp$Blobs')
pylab.title('Number of vortex blobs')