"""

    ex_03
    
        Example 03 showing the usage of plotting routines
        
        Generate (nBlobs x nBlobs) vortex blobs evenly spaced in the interval 
        [-1, 1] x [-1, 1] with associated circulation given by:
            
            w = 100.0* abs(deltaX*deltaX*sin(x))
        
        with:
            
            deltaX = deltaY = 2.0/(nBlobs-1.0)
        
        Plot the vorticity using circulation and induced vorticity and pyplot
        and dolfin plotting systems.

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
# First added:  2013-07-01                                                                                                          
# Last changed: 2013-07-01
# -*- coding: utf-8 -*-

# -*- coding: utf-8 -*-

"""    
    Reviews:
        
"""

# required packages
import numpy
import pylab
import pHyFlow
from mpl_toolkits.mplot3d import Axes3D

#=============================================================================
# input parameters

# parameters definition
nBlobs = 55         # the number of vortex blobs
nTargets = 22       # the number of target points where to evaluate velocities

overlap = 0.5       # the overlap between vortex blobs (together with the number
                    # of particles defines the core size of the particles)

k = 2.0;            # the peakedness of the kernel typical value is 2.0

blockSize = 64      # the size of the block of the GPU, this value is strongly
                    # related to the speed of the calculations on the GPU, for
                    # a large number of particles (> 1024) 128 is a good number.

xBounds = numpy.array([-1,3]) # define the bounds for the plot
yBounds = numpy.array([-1,1]) # define the bounds for the plot

nPlotPoints = [200,100]

def wExact(x,y): return (10.0/(4.0*numpy.pi*0.01*1.0))*numpy.exp(-(x*x+y*y)/(4.0*0.01*1.0))
#=============================================================================



#=============================================================================
# computation of parameters

h = 2.0/(nBlobs-1.0)   # spacing in x of vortex blobs (note that particles
                            # have the same spacing in y => deltaY = deltaX)

sigma = (h/overlap) # compute the core spreading
sigmasqr = sigma * sigma # compute the square of the core spreading
#=============================================================================



#=============================================================================
# computation of vortex blobs

# generate the coordinates of the vortex blobs
xB = numpy.linspace(-1.0,1.0,nBlobs)        # nBlobs evenly spaced in x \in [-1, 1]
yB = numpy.linspace(-1.0,1.0,nBlobs)        # nBlobs evenly spaced in y \in [-1, 1]
xB,yB = numpy.meshgrid(xB,yB)               # generate the 2d grid of blobs in [-1, 1] x [-1, 1]

# flatten the grids into vectors
xB = xB.flatten()                           
yB = yB.flatten()

# compute the circulations of the blobs
wB = wExact(xB,yB)*h*h
#=============================================================================



#=============================================================================
# plot vorticity from circulation using pyplot
figureHandle = pHyFlow.vortex.plotvorticity(xBounds,yBounds,xB,yB,wB,\
                                            sigma,overlap,k=2,\
                                            wType=pHyFlow.options.WBLOB_VORTICITY,\
                                            figureHandle=None,plotType=pHyFlow.options.PYLAB_PLOT)
                                            
# plot vorticity from circulation using dolfin
figureHandle = pHyFlow.vortex.plotvorticity(xBounds,yBounds,xB,yB,wB,\
                                            sigma,overlap,k=2,\
                                            wType=pHyFlow.options.WBLOB_VORTICITY,\
                                            figureHandle=None,plotType=pHyFlow.options.DOLFIN_PLOT,interactive=True)

# plot vorticity from vorticity function using pyplot
figureHandle = pHyFlow.vortex.plotvorticity(xBounds,yBounds,xB,yB,wB,sigma,overlap,k=2,kernel=pHyFlow.options.GAUSS_KERNEL,\
                                            nPlotPoints=nPlotPoints,hardware=pHyFlow.options.GPU_HARDWARE,\
                                            blocksize=128,method=pHyFlow.options.DIRECT_METHOD,\
                                            wType=pHyFlow.options.WFUNCTION_VORTICITY,figureHandle=None,\
                                            plotType=pHyFlow.options.PYLAB_PLOT)
                                            
# plot vorticity from vorticity function using dolfin
figureHandle = pHyFlow.vortex.plotvorticity(xBounds,yBounds,xB,yB,wB,\
                                            sigma,overlap,k=2,kernel=pHyFlow.options.GAUSS_KERNEL,\
                                            nPlotPoints=nPlotPoints,hardware=pHyFlow.options.GPU_HARDWARE,\
                                            blocksize=128,method=pHyFlow.options.DIRECT_METHOD,\
                                            wType=pHyFlow.options.WFUNCTION_VORTICITY,\
                                            figureHandle=None,plotType=pHyFlow.options.DOLFIN_PLOT,interactive=True)                                            
                                            
# plot the exact vorticity (using imshow)
pylab.figure()
pylab.imshow(wExact(xB.reshape([nBlobs,nBlobs]),yB.reshape([nBlobs,nBlobs])),extent=(xB.min()-0.5*sigma*overlap,\
                           xB.max()+0.5*sigma*overlap,\
                           yB.min()-0.5*sigma*overlap,\
                           yB.max()+0.5*sigma*overlap),\
                           interpolation='nearest',origin='lower')
                           
figure3d = pylab.figure()
ax = Axes3D(figure3d)
ax.plot_surface(xB.reshape([nBlobs,nBlobs]),yB.reshape([nBlobs,nBlobs]),wExact(xB.reshape([nBlobs,nBlobs]),yB.reshape([nBlobs,nBlobs])),rstride=1, cstride=1)

# plot the exact vorticity (using pcolormesh)
# generate the coordinates of the plot points
# this is done because pcolor uses the corner points of the input coordinates
# to define the rectangle, not the center point
xP = numpy.linspace(-1.0-h,1.0+h,nBlobs+1)        # nBlobs+1 points evenly spaced in x \in [-1, 1]
yP = numpy.linspace(-1.0-h,1.0+h,nBlobs+1)        # nBlobs+1 points evenly spaced in y \in [-1, 1]
xP,yP = numpy.meshgrid(xP,yP)               # generate the 2d grid of plot points in [-1, 1] x [-1, 1]                           

pylab.figure()
pylab.pcolormesh(xP,yP,wExact(xB.reshape([nBlobs,nBlobs]),yB.reshape([nBlobs,nBlobs])))
#=============================================================================