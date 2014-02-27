"""

Test the navier-stokes modules: _navierStokes.py

"""

import pHyFlow

# External modules
import dolfin
import numpy as np
import pylab as py
py.ion()


# Geometry Position
cmGlobal = np.array([0.,0.])
thetaLocal = np.deg2rad(0.)

# Fluid Parameters
Re = 5000.
uInf = 1.0
uMax = 6.0
L   = 1.0
nu   = (uInf*L)/(Re) # Unit Geometry with Re = 1000.

# CFL Number
cfl = 0.95

# Probe Grid Parameters
x0,y0,Lx,Ly,nx,ny = -1.0,-1.0,2.0,2.0,10,10    
#-----------------------------------------------------------------------------




#-----------------------------------------------------------------------------   
# Mmesh and boundary mesh files

# Export to files
meshFileName = './geometry/fullGrid_NACA0012_nVertices-57k_pressureOutlet.xml.gz'
boundaryDomainsFileName = './geometry/fullGrid_NACA0012_nVertices-57k_pressureOutlet_facet_region.xml.gz'

#-----------------------------------------------------------------------------   


#------------------------------------------------------------------------------
# Cylinder 2D test case

geometry = {'mesh' : meshFileName,
            'boundaryDomains': boundaryDomainsFileName,
            'cmGlobal': cmGlobal,
            'thetaLocal': thetaLocal}

probeGrid = {'origin': np.array([x0,y0]),
             'L'     : np.array([Lx,Ly]),
             'N'     : np.array([nx,ny])}            

# Initialize Navier-Stokes problem
NSDomain = pHyFlow.navierStokes.NavierStokes(geometry,probeGrid,uMax,nu,cfl)
                                     
#------------------------------------------------------------------------------
                                   
                                   
#------------------------------------------------------------------------------                                     
                                     
# Get boundary coordinates
x,y = NSDomain.getCoordinates()
py.figure()
py.plot(x,y,'b.')
py.title('Navier-Stokes Vector DOFs')
py.xlabel(r'$x$')
py.ylabel(r'$y$')
py.show(block=True)
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Get dirichlet boundary coordinates
xB,yB = NSDomain.getBoundaryCoordinates()
py.figure()
py.plot(x,y,'b.')
py.plot(xB,yB,'ko')
py.title('Navier-Stokes Dirichlet boundary Vector DOFs')
py.xlabel(r'$x$')
py.ylabel(r'$y$')
py.show(block=True)
#------------------------------------------------------------------------------

# Evolve NSDomain
vxBoundary, vyBoundary = 1.0*np.ones(xB.shape), np.zeros(yB.shape)
cmGlobalNew, thetaGlobalNew = cmGlobal,thetaLocal
cmDotGlobal, thetaGlobal = np.array([0.,0.]), 0.

saveDir = './data/'
vFile = dolfin.File(saveDir + "velocity.pvd", "compressed")
pFile = dolfin.File(saveDir + "pressure.pvd", "compressed")
wFile = dolfin.File(saveDir + "vorticity.pvd", "compressed")

solver = NSDomain._NavierStokes__solver

for i in xrange(int(1e6)):

    T = i*NSDomain.deltaT

    NSDomain.evolve(vxBoundary,vyBoundary,cmGlobalNew,thetaGlobalNew,
                    cmDotGlobal,thetaGlobal)
    
    diff_u = dolfin.norm(solver.u1) - dolfin.norm(solver.u0)            

    if i % 5 == 0:                
        #dolfin.plot(dolfin.sqrt(dolfin.inner(NSDomain._solver.u1,NSDomain._solver.u1)),key='vNorm')         
        vFile << (NSDomain._NavierStokes__solver.u1, T)
        pFile << (NSDomain._NavierStokes__solver.p1, T)
        wFile << (NSDomain._NavierStokes__solver.vorticity(), T)

    print "Step = %g\tT = %g" % (i, T)
    print "Difference in velocity : %g" % diff_u    
    