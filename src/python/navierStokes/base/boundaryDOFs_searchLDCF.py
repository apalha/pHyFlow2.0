# -*- coding: utf-8 -*-
"""
Created on Thu Feb 13 14:41:29 2014

@author: lento
"""

import dolfin
import numpy
import time
from pHyFlow.navierStokes.base import boundary


#-----------------------------------------------------------------------------    
# Lid Driven Cavity Flow - Boundary Conditions

N = 200

# Make mesh
mesh = dolfin.UnitSquareMesh(N,N)

V = dolfin.VectorFunctionSpace(mesh,'CG',2)
    
 # Sub domain for no-slip (mark whole boundary, Lid will overwrite)
class NoSlip(dolfin.SubDomain):
    def inside(self,x,on_boundary):
        return on_boundary

# Sub domain for Lid (top)
class Lid(dolfin.SubDomain):
    def inside(self,x,on_boundary):
        return x[1] > 1.0 - dolfin.DOLFIN_EPS and on_boundary
  
# Define boundary domains    
boundaryDomains = dolfin.MeshFunction('size_t', mesh, mesh.topology().dim()-1)

boundaryDomains.set_all(1)

# Mark no-slip
noslip = NoSlip()
noslip.mark(boundaryDomains, 2) # No-slip ID: 2

# Mark exterior I.D
lid = Lid()
lid.mark(boundaryDomains, 3) # extID : 3  
  
#-----------------------------------------------------------------------------   
  
#-----------------------------------------------------------------------------   
# SLOW, FOR loop search

# Determine both 
start = time.time()
boundary_DOFCoordinates, boundary_vectorDOFIndex = boundary.locate_boundaryDOFs_slow(mesh,boundaryDomains,3)
S1 = (time.time()-start)
print "Slow evaluate (Total) :: %g " % S1
#-----------------------------------------------------------------------------     
  
#-----------------------------------------------------------------------------   
# FAST, SIMPLER method
# 

# Determine the vector boundary DOF indicies
start = time.time()
xyIndices = boundary.vectorDOF_boundaryIndex(V,boundaryDomains,3)
F1 = (time.time()-start)
print "Fast evaluate (find vector DOF boundary indices) :: %g " % F1

# Determine the vector boundary DOF indicies
start = time.time()
xyBoundaryCoordinates = boundary.vectorDOF_coordinates(mesh,V,xyIndices[0])
F2 = (time.time()-start)
print "Fast evaluate (find vector DOF boundary coordinates) :: %g " % F2
print "Fast evaluate (Total) :: %g " % (F1+F2)

print "Speed gain: %gx with %gx%g cells." % (S1/(F1+F2),N,N)



#-----------------------------------------------------------------------------   
