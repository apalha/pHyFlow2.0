"""

    ex_04
    
        Example 04 showing the usage of time evolution routines
        
        Generate two vortex blobs at V1 = (0.0,-1.0) and V2 = (0.0, 1.0)
        with circulations G1 = -50.0 and G2 = 50.0. Total vorticity field
        given by:
            
            w(x,y) = G1 * exp(-((x^2)+(y-1.0)^2)/(2*sigmasqr))/(2*pi*sigmasqr) +
                     + G2 * exp(-((x^2)+(y+1.0)^2)/(2*sigmasqr))/(2*pi*sigmasqr)
            
        with:
            
            sigmasqr the square of the core radius
        
        Advancing in time the vortex blobs without remeshing. Compare evolution
        with forward Euler and RK4.

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
# First added:  2013-07-19                                                                                                          

"""    
    Reviews:
        
"""

import pHyFlow
import numpy
from matplotlib import rc
import pylab

pylab.ion() # turn interactive mode on to show animation

rc('text',usetex=True) # to use TeX in pyplot

#=============================================================================
# input parameters

# parameters definition
sigma = 0.01        # the core radius of the vortex blobs

overlap = 0.5       # the overlap between vortex blobs

k = 2.0;            # the peakedness of the kernel typical value is 2.0

blockSize = 128      # the size of the block of the GPU, this value is strongly
                    # related to the speed of the calculations on the GPU, for
                    # a large number of particles (> 1024) 128 is a good number.

dt = 0.01           # the time step
tSteps = 200        # the number of time steps

# generate the coordinates of the vortex blobs
xB = numpy.array([0.0,0.0])
yB = numpy.array([1.0,-1.0]) 

# generate the vorticity of the blobs
wB = numpy.array([50.0,50.0])

# define the domains for plotting and for computation
xBoundsDomain = numpy.array([-2.0,2.0]) # define the bounds for the domain
yBoundsDomain = numpy.array([-2.0,2.0]) # define the bounds for the domain

xBoundsPlot = numpy.array([-2.0,2.0]) # define the bounds for the plot
yBoundsPlot = numpy.array([-2.0,2.0]) # define the bounds for the plot

nPlotPoints = [100,100] # the number of plot points

#=============================================================================


#=============================================================================
# computation of parameters

sigmasqr = sigma * sigma # compute the square of the core spreading
h = sigma*overlap # compute the spacing of the blobs
#=============================================================================


#=============================================================================
# set initial positions for forward Euler and Runge-Kutta 4th order integrators
xB_FE = numpy.zeros([tSteps+1,xB.size])
yB_FE = numpy.zeros([tSteps+1,yB.size])
wB_FE = numpy.zeros([tSteps+1,wB.size])

xB_FE[0,:] = xB.copy()
yB_FE[0,:] = yB.copy()
wB_FE[0,:] = wB.copy()

xB_RK4 = numpy.zeros([tSteps+1,xB.size])
yB_RK4 = numpy.zeros([tSteps+1,yB.size])
wB_RK4 = numpy.zeros([tSteps+1,wB.size])

xB_RK4[0,:] = xB.copy()
yB_RK4[0,:] = yB.copy()
wB_RK4[0,:] = wB.copy()

#=============================================================================


#=============================================================================
# plot vorticity and blobs
                                            
# plot vorticity from vorticity function using dolfin (forward Euler)
figureHandleBlobs_FE = pHyFlow.vortex.plotvorticity(xBoundsPlot,yBoundsPlot,xB_FE[0,:],yB_FE[0,:],wB_FE[0,:],\
                                            sigma,overlap,k=k,kernel=pHyFlow.options.GAUSS_KERNEL,\
                                            nPlotPoints=nPlotPoints,hardware=pHyFlow.options.GPU_HARDWARE,\
                                            blocksize=128,method=pHyFlow.options.DIRECT_METHOD,\
                                            wType=pHyFlow.options.WFUNCTION_VORTICITY,\
                                            figureHandle=None,title='Vorticity from blob function (FE)',plotType=pHyFlow.options.DOLFIN_PLOT,interactive=True) 

# plot vorticity from vorticity function using dolfin (Runge-Kutta 4th order)
figureHandleBlobs_RK4 = pHyFlow.vortex.plotvorticity(xBoundsPlot,yBoundsPlot,xB_RK4[0,:],yB_RK4[0,:],wB_RK4[0,:],\
                                            sigma,overlap,k=k,kernel=pHyFlow.options.GAUSS_KERNEL,\
                                            nPlotPoints=nPlotPoints,hardware=pHyFlow.options.GPU_HARDWARE,\
                                            blocksize=128,method=pHyFlow.options.DIRECT_METHOD,\
                                            wType=pHyFlow.options.WFUNCTION_VORTICITY,\
                                            figureHandle=None,title='Vorticity from blob function (RK4)',plotType=pHyFlow.options.DOLFIN_PLOT,interactive=True) 

# plot the vortex blob centers
pylab.plot(xB_RK4[0,0],yB_RK4[0,0],'-r',xB_RK4[0,1],yB_RK4[0,1],'-r',xB_RK4[0,:],yB_RK4[0,:],'or',\
           xB_FE[0,0],yB_FE[0,0],'-b',xB_FE[0,0],xB_FE[0,1],yB_FE[0,1],'-b',xB_FE[0,:],yB_FE[0,:],'or')

pylab.axis('equal')
pylab.axis(numpy.concatenate((xBoundsPlot,yBoundsPlot)))

pylab.draw()

#=============================================================================


#=============================================================================
# advance blobs in time

# allocate memory space for the core separation radius
r_RK4 = numpy.zeros(tSteps+1) # Runge-Kutta 4th order radius
r_FE = numpy.zeros(tSteps+1) # forward Euler radius

r_RK4[0] = numpy.sqrt(((xB[1]-xB[0])**2)+((yB[1]-yB[0])**2)) # compute the initial radius between the two blobs
r_FE[0] = r_RK4[0]

for tStep in range(1,tSteps+1):
    print '----------------------------------------------------------'
    print '   Time step :: ' + str(tStep)

    # move vortex blobs (FE)   
    xB_FE[tStep,:],yB_FE[tStep,:],wB_FE[tStep,:] = pHyFlow.vortex.evolution.convection(dt,xB_FE[tStep-1,:],yB_FE[tStep-1,:],wB_FE[tStep-1,:],sigma,k=k,kernel=pHyFlow.options.GAUSS_KERNEL,\
                                                   hardware=pHyFlow.options.GPU_HARDWARE,method=pHyFlow.options.DIRECT_METHOD,\
                                                   blocksize=128,integrator=pHyFlow.options.FE_INTEGRATOR)
                                                   
    # move vortex blobs (RK4)   
    xB_RK4[tStep,:],yB_RK4[tStep,:],wB_RK4[tStep,:] = pHyFlow.vortex.evolution.convection(dt,xB_RK4[tStep-1,:],yB_RK4[tStep-1,:],wB_RK4[tStep-1,:],sigma,k=k,kernel=pHyFlow.options.GAUSS_KERNEL,\
                                                   hardware=pHyFlow.options.GPU_HARDWARE,method=pHyFlow.options.DIRECT_METHOD,\
                                                   blocksize=128,integrator=pHyFlow.options.RK4_INTEGRATOR)                                                   
    
    # plot vorticity distribution (FE)
    figureHandleBlobs_FE = pHyFlow.vortex.plotvorticity(xBoundsPlot,yBoundsPlot,xB_FE[tStep,:],yB_FE[tStep,:],wB_FE[tStep,:],\
                                            sigma,overlap,k=k,\
                                            wType=pHyFlow.options.WFUNCTION_VORTICITY,\
                                            figureHandle=figureHandleBlobs_FE,title='Vorticity from blob function (FE)',\
                                            nPlotPoints=nPlotPoints,plotType=pHyFlow.options.DOLFIN_PLOT,interactive=True)      

    # plot vorticity distribution (RK4)
    figureHandleBlobs_RK4 = pHyFlow.vortex.plotvorticity(xBoundsPlot,yBoundsPlot,xB_RK4[tStep,:],yB_RK4[tStep,:],wB_RK4[tStep,:],\
                                            sigma,overlap,k=k,\
                                            wType=pHyFlow.options.WFUNCTION_VORTICITY,\
                                            figureHandle=figureHandleBlobs_RK4,title='Vorticity from blob function (RK4)',\
                                            nPlotPoints=nPlotPoints,plotType=pHyFlow.options.DOLFIN_PLOT,interactive=True)                                               
    
    # plot the vortex blob centers
    pylab.cla()
    pylab.plot(xB_RK4[0:tStep+1,0],yB_RK4[0:tStep+1,0],'-r',\
               xB_FE[0:tStep+1,0],yB_FE[0:tStep+1,0],'-b',\
               xB_RK4[0:tStep+1,1],yB_RK4[0:tStep+1,1],'-r',\
               xB_RK4[tStep,:],yB_RK4[tStep,:],'or',\
               xB_FE[0:tStep+1,1],yB_FE[0:tStep+1,1],'-b',\
               xB_FE[tStep,:],yB_FE[tStep,:],'ob')

    pylab.axis('equal')    
    pylab.axis(numpy.concatenate((xBoundsPlot,yBoundsPlot)))
    pylab.legend(['RK4','FE'],loc='upper right',shadow=True)                        
    pylab.draw()                           

    # compute the radius between the two blobs                                            
    r_FE[tStep]= numpy.sqrt(((xB_FE[tStep,1]-xB_FE[tStep,0])**2)+((yB_FE[tStep,1]-yB_FE[tStep,0])**2)) # forward Euler
    r_RK4[tStep]= numpy.sqrt(((xB_RK4[tStep,1]-xB_RK4[tStep,0])**2)+((yB_RK4[tStep,1]-yB_RK4[tStep,0])**2)) # RK4
    
    print '                                   FE               RK4'
    print '   Core separation       :: ' + '%14.10f  %14.10f' % (r_FE[tStep], r_RK4[tStep])
    print '   Core separation error :: ' + '%14.10f  %14.10f' % (r_FE[tStep]-r_FE[0], r_RK4[tStep]-r_RK4[0])
    print '----------------------------------------------------------'
    
#=============================================================================

# plot the error evolution in time
pylab.figure()
pylab.semilogy(numpy.arange(0,tSteps+1)*dt,numpy.abs(r_FE-r_FE[0]),'o-b',numpy.arange(0,tSteps+1)*dt,numpy.abs(r_RK4-r_RK4[0]),'o-r')
pylab.legend(['FE', 'RK4'],loc='lower right',shadow = True)
pylab.xlabel(r't [s]')
pylab.ylabel(r'$|r(t)-r(0)|$')
pylab.draw()