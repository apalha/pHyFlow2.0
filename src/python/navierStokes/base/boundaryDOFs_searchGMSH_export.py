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

nVertices = '13k'

# Make mesh
mesh = dolfin.Mesh('cylinder2D_nVertices-13k_delQuad_noBL_finerBody.xml.gz')

V = dolfin.VectorFunctionSpace(mesh,'CG',2)
  
# Define boundary domains    
boundaryDomains = dolfin.MeshFunction('size_t',mesh,'cylinder2D_nVertices-13k_delQuad_noBL_finerBody_facet_region.xml.gz')
  
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

print "Speed gain: %gx with %s vertices." % (S1/(F1+F2),nVertices)



#-----------------------------------------------------------------------------   
