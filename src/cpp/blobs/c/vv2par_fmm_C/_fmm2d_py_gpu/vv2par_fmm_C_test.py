# -*- coding: utf-8 -*-
"""
Created on Mon Sep 23 14:29:09 2013

@author: gorkiana
"""

import numpy
import pyublas
from _fmm2d_py import fmm2d_py
import pylab
import time
import fmm_velocity

pylab.ion()
#------------------------------------------------------------------------------
# Generate exact vorticity and blobs in a regular grid

# the vorticity
w_exact = lambda x,y,gamma,pVortex,epsilon: gamma*(numpy.exp(-(((x-pVortex[0])**2)\
                                            +((y-pVortex[1])**2))/(2.0*epsilon*epsilon)))/(2.0*numpy.pi*epsilon*epsilon)

# the parameters of the vorticity
gamma = 50.0
epsilon = 0.5*numpy.sqrt(numpy.pi/12.0)
pVortex = [0.0, 0.0]

# generate the blobs evenly spaced
nBlobs = 200
overlap = 0.5

x = numpy.linspace(-3.0,3.0,nBlobs)
y = x.copy()

h = x[1]-x[0] # compute the size of the blobs

sigma = h/overlap # compute the size of the blobs

[xBlob,yBlob] = numpy.meshgrid(x,y) # compute the coordinates of all the blobs in a regular grid

xBlob = xBlob.T.flatten().copy() # flatten them out (the .T for transpose is used because we wish to compare with Matlab results)
yBlob = yBlob.T.flatten().copy()

wBlob = w_exact(xBlob,yBlob,gamma,pVortex,epsilon)*h*h

#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Plot the blobs

pylab.scatter(xBlob,yBlob,c=wBlob,s=500,linewidths=0.0,alpha=0.7)

#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Compute the induced velocities

# convert blob parameters into the parameters of the fmm solver
alpha = 1.2564312086261696770
k = 2.0
ksigmasqr = k*sigma*sigma
xopt = numpy.sqrt(alpha*ksigmasqr)
cutoff = 5.0*xopt;
Ndirect = 35
tol = 1.0e-6

xTarget = xBlob.copy()
yTarget = yBlob.copy()

start = time.time()
vx_fmm_cpu,vy_fmm_cpu = fmm2d_py(xBlob,yBlob,wBlob,xTarget,yTarget,Ndirect,xopt,cutoff,tol)
print "time to compute fmm computation :: " + str(time.time()-start)

start = time.time()
vx_direct_cpu,vy_direct_cpu = fmm2d_py(xBlob,yBlob,wBlob,xTarget,yTarget,Ndirect,xopt,cutoff,0.0)
print "time to compute direct computation :: " + str(time.time()-start)

#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Plot the induced velocities
pylab.figure()
pylab.quiver(xTarget,yTarget,-vy_fmm_cpu,-vx_fmm_cpu) # the fmm velocity
pylab.quiver(xTarget,yTarget,-vy_direct_cpu,-vx_direct_cpu) # the direct velocity

pylab.figure()
pylab.quiver(xTarget,yTarget,-vy_fmm_cpu+vy_direct_cpu,-vx_fmm_cpu+vx_direct_cpu,scale=1.0) # the difference between fmm and direct

#------------------------------------------------------------------------------

vx,vy = fmm_velocity.fmm_velocity(xBlob,yBlob,wBlob,sigma,k=2.0,xTarget=xTarget,yTarget=yTarget,Ndirect=Ndirect,tol=1.0e-6)












