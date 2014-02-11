"""

    ex_02
    
        Example 02 showing the usage of functions to compute induced vorticity
        using different kernels and hardware (CPU and GPU).
        
        Generate (nBlobs x nBlobs) vortex blobs evenly spaced in the interval 
        [-1, 1] x [-1, 1] with associated circulation given by:
            
            w = 100.0* abs(deltaX*deltaX*sin(x))
        
        with:
            
            deltaX = deltaY = 2.0/(nBlobs-1.0)
        
        Generate (nTargets x nTargets) target points evenly spaced in the
        interval [-1, 1] x [-1, 1], where to compute the induced velocities.
        
        Compare times for computing with CPU and GPU using direct calculation.

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
# First added:  2013-06-26                                                                                                          
# Last changed: 2013-06-26
# -*- coding: utf-8 -*-

# -*- coding: utf-8 -*-

"""    
    Reviews:
        
"""

# optional packages used just to determine the computation times of each routine
# for comparison purposes
import time

# required packages
import numpy
import pHyFlow

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
#=============================================================================



#=============================================================================
# computation of parameters

deltaX = 2.0/(nBlobs-1.0)   # spacing in x of vortex blobs (note that particles
                            # have the same spacing in y => deltaY = deltaX)

sigmasqr = (deltaX/overlap)*(deltaX/overlap) # compute the square of the core spreading
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
wB = 100.0*numpy.abs(deltaX*deltaX*numpy.sin(xB))
#=============================================================================



#=============================================================================
# computation of target points

# compute the coordinates of the target points
xT = numpy.linspace(-1.0,1.0,nTargets)  # nTargets evenly spaced in x \in [-1, 1]
yT = numpy.linspace(-1.0,1.0,nTargets)  # nTargets evenly spaced in y \in [-1, 1]
xT,yT = numpy.meshgrid(xT,yT)           # generate the 2d grid of target points in [-1, 1] x [-1, 1]

# flatten the grids into vectors
xT = xT.flatten()                           
yT = yT.flatten()
#=============================================================================


#=============================================================================
# computation of induced vorticity

#----------------------------------------------------------------------------
# Gauss kernel GPU, direct calculation

# run computation of induced vorticity once without checking time (the first time is frequently much slower)
wTemp = pHyFlow.vortex.induced.vorticity(xB,yB,wB,numpy.sqrt(sigmasqr),k=k,kernel=pHyFlow.options.GAUSS_KERNEL,\
                                         xEval=xT,yEval=yT,hardware=pHyFlow.options.GPU_HARDWARE,\
                                         blocksize=blockSize,method=pHyFlow.options.DIRECT_METHOD)

# start time counter
startTime_ggd = time.time() # ggd: Gauss Gpu Direct

# run computation of induced vorticity: Gauss kernel, GPU and direct calculation
w_ggd = pHyFlow.vortex.induced.vorticity(xB,yB,wB,numpy.sqrt(sigmasqr),k=k,kernel=pHyFlow.options.GAUSS_KERNEL,\
                                         xEval=xT,yEval=yT,hardware=pHyFlow.options.GPU_HARDWARE,\
                                         blocksize=blockSize,method=pHyFlow.options.DIRECT_METHOD)

# stop time counter                                                
endTime_ggd = time.time()                                                

# display time
print 'Time it took for Gauss kernel, GPU and direct ::' + str(endTime_ggd-startTime_ggd)

#----------------------------------------------------------------------------


#----------------------------------------------------------------------------
# Gauss kernel CPU, direct calculation

# run computation of induced vorticity once without checking time (the first time is frequently much slower)
wTemp = pHyFlow.vortex.induced.vorticity(xB,yB,wB,numpy.sqrt(sigmasqr),k=k,kernel=pHyFlow.options.GAUSS_KERNEL,\
                                         xEval=xT,yEval=yT,hardware=pHyFlow.options.CPU_HARDWARE,\
                                         blocksize=blockSize,method=pHyFlow.options.DIRECT_METHOD)

# start time counter
startTime_gcd = time.time() # gcd: Gauss Cpu Direct

# run computation of induced vorticity: Gauss kernel, CPU and direct calculation
w_gcd = pHyFlow.vortex.induced.vorticity(xB,yB,wB,numpy.sqrt(sigmasqr),k=k,kernel=pHyFlow.options.GAUSS_KERNEL,\
                                         xEval=xT,yEval=yT,hardware=pHyFlow.options.CPU_HARDWARE,\
                                         blocksize=blockSize,method=pHyFlow.options.DIRECT_METHOD)

# stop time counter                                                
endTime_gcd = time.time()                                                

# display time
print 'Time it took for Gauss kernel, CPU and direct ::' + str(endTime_gcd-startTime_gcd)

#----------------------------------------------------------------------------



#----------------------------------------------------------------------------
# Cutoff kernel GPU, direct calculation

# run computation of induced vorticity once without checking time (the first time is frequently much slower)
wTemp = pHyFlow.vortex.induced.vorticity(xB,yB,wB,numpy.sqrt(sigmasqr),k=k,kernel=pHyFlow.options.CUTOFF_KERNEL,\
                                         xEval=xT,yEval=yT,hardware=pHyFlow.options.GPU_HARDWARE,\
                                         blocksize=blockSize,method=pHyFlow.options.DIRECT_METHOD)

# start time counter
startTime_cgd = time.time() # cgd: Cutoff Gpu Direct

# run computation of induced vorticity: Gauss kernel, GPU and direct calculation
w_cgd = pHyFlow.vortex.induced.vorticity(xB,yB,wB,numpy.sqrt(sigmasqr),k=k,kernel=pHyFlow.options.CUTOFF_KERNEL,\
                                         xEval=xT,yEval=yT,hardware=pHyFlow.options.GPU_HARDWARE,\
                                         blocksize=blockSize,method=pHyFlow.options.DIRECT_METHOD)

# stop time counter                                                
endTime_cgd = time.time()                                                

# display time
print 'Time it took for Cutoff kernel, GPU and direct ::' + str(endTime_cgd-startTime_cgd)

#----------------------------------------------------------------------------


#----------------------------------------------------------------------------
# Cutoff kernel CPU, direct calculation

# run computation of induced vorticity once without checking time (the first time is frequently much slower)
wTemp = pHyFlow.vortex.induced.vorticity(xB,yB,wB,numpy.sqrt(sigmasqr),k=k,kernel=pHyFlow.options.CUTOFF_KERNEL,\
                                         xEval=xT,yEval=yT,hardware=pHyFlow.options.CPU_HARDWARE,\
                                         blocksize=blockSize,method=pHyFlow.options.DIRECT_METHOD)

# start time counter
startTime_ccd = time.time() # ccd: Cutoff Cpu Direct

# run computation of induced velocities: Cutoff kernel, CPU and direct calculation
w_ccd = pHyFlow.vortex.induced.vorticity(xB,yB,wB,numpy.sqrt(sigmasqr),k=k,kernel=pHyFlow.options.CUTOFF_KERNEL,\
                                         xEval=xT,yEval=yT,hardware=pHyFlow.options.CPU_HARDWARE,\
                                         blocksize=blockSize,method=pHyFlow.options.DIRECT_METHOD)

# stop time counter                                                
endTime_ccd = time.time()                                                

# display time
print 'Time it took for Cutoff kernel, CPU and direct ::' + str(endTime_ccd-startTime_ccd)

#----------------------------------------------------------------------------
#=============================================================================



#=============================================================================
# computation of errors between GPU and CPU

print '\n'

print 'Max error in the two functions (Gaussian for GPU vs CPU) :: \n' +\
      '    e_w = ' + str(numpy.abs(w_ggd - w_gcd).max())

print '\n'

print 'Max error in the two functions (CUTOFF for GPU vs CPU) :: \n' +\
      '    e_w = ' + str(numpy.abs(w_cgd - w_ccd).max())
#=============================================================================