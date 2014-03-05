"""

Test the navier-stokes modules: _navierStokes.py

"""

import pHyFlow

# External modules
import dolfin
import numpy as np
import pylab as py
py.ion()


#-----------------------------------------------------------------------------    

# Mesh Data
N = 32 
#mesh, boundaryDomains = lidDrivenCavityFlow(N)

# Geometry Position
cmGlobal = np.array([0.,0.])
thetaLocal = np.deg2rad(0.)

# Fluid Parameters
Re = 1.
uMax = 1.0
nu   = (uMax*1.0)/(Re) # Unit Geometry with Re = 1000.

# CFL Number
cfl = 0.95

# Probe Grid Parameters
x0,y0,Lx,Ly,nx,ny = 0.2,0.2,0.6,0.6,10,10    
#-----------------------------------------------------------------------------


#-----------------------------------------------------------------------------
# Define Mesh
mesh = dolfin.UnitSquareMesh(N,N)

#Non-uniform mesh
xy = mesh.coordinates()
xyi = 0.5-0.5*np.cos(np.pi*xy)
mesh.coordinates()[:] = xyi

#-----------------------------------------------------------------------------
    
#-----------------------------------------------------------------------------    
# Lid Driven Cavity Flow - Boundary Conditions
    
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
noslip.mark(boundaryDomains, pHyFlow.eulerian.eulerianOptions.ID_NOSLIP_BOUNDARY) # No-slip ID: 2

# Mark exterior I.D
lid = Lid()
lid.mark(boundaryDomains, pHyFlow.eulerian.eulerianOptions.ID_EXTERNAL_BOUNDARY) # extID : 3  
  
#-----------------------------------------------------------------------------   


#-----------------------------------------------------------------------------   
# Export mesh and boundary mesh files

# Export to files
meshFileName = 'LDCF_%gx%g_mesh.xml.gz' % (N,N)
boundaryDomainsFileName = 'LDCF_%gx%g_boundaryDomains.xml.gz' % (N,N)
dolfin.File(meshFileName) << mesh
dolfin.File(boundaryDomainsFileName) << boundaryDomains
#-----------------------------------------------------------------------------   


#------------------------------------------------------------------------------
# Lid-Driven cavity flow test case

geometry = {'mesh' : meshFileName,
            'boundaryDomains': boundaryDomainsFileName,
            'cmGlobal': cmGlobal,
            'thetaLocal': thetaLocal}

probeGrid = {'origin': np.array([x0,y0]),
             'L'     : np.array([Lx,Ly]),
             'N'     : np.array([nx,ny])}            

# Initialize Navier-Stokes problem
eulerian = pHyFlow.eulerian.EulerianSolver(geometry,probeGrid,
                                           uMax=uMax,nu=nu,cfl=cfl,
                                           deltaT=1e-5)
                                     
#------------------------------------------------------------------------------
                                     
#------------------------------------------------------------------------------                                     
                                     
# Get boundary coordinates
x,y = eulerian.getCoordinates()
py.figure()
py.plot(x,y,'b.')
py.title('Navier-Stokes Vector DOFs')
py.xlabel(r'$x$')
py.ylabel(r'$y$')
py.show(block=True)
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Plot Velocity vectors - ONLY FOR TESTING (else need to output to pvd)
vx,vy = np.random.rand(2,x.shape[0])
eulerian.setVelocity(vx,vy)
dolfin.plot(eulerian._EulerianSolver__solver.u1, title='ONLY FOR TESTING - Random velocity field')
dolfin.interactive()

#------------------------------------------------------------------------------
# Set velocity
vx,vy = np.zeros((2,x.shape[0]))
eulerian.setVelocity(vx,vy)
# PLlot ONLY FOR TESTING (else need to output to pvd)
dolfin.plot(eulerian._EulerianSolver__solver.u1, title='ONLY FOR TESTING - Changed to zero velocity field')
dolfin.interactive()
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Get dirichlet boundary coordinates
xB,yB = eulerian.getBoundaryCoordinates()
py.figure()
py.plot(x,y,'b.')
py.plot(xB,yB,'ko')
py.title('Navier-Stokes Dirichlet boundary Vector DOFs')
py.xlabel(r'$x$')
py.ylabel(r'$y$')
py.show(block=True)
#------------------------------------------------------------------------------

# Evolve NSDomain
vxBoundary, vyBoundary = -1.0*np.ones(xB.shape), np.zeros(yB.shape)
cmGlobalNew, thetaGlobalNew = cmGlobal,thetaLocal
cmDotGlobal, thetaGlobal = np.array([0.,0.]), 0.

saveDir = './data/'
vFile = dolfin.File(saveDir + "velocity.pvd", "compressed")
pFile = dolfin.File(saveDir + "pressure.pvd", "compressed")
wFile = dolfin.File(saveDir + "vorticity.pvd", "compressed")

solver = eulerian._EulerianSolver__solver

for timeStep in range(1,100):

    T = timeStep*eulerian.deltaT

    eulerian.evolve(vxBoundary,vyBoundary,cmGlobalNew,thetaGlobalNew,
                    cmDotGlobal,thetaGlobal)
    
    diff_u = dolfin.norm(solver.u1) - dolfin.norm(solver.u0)            

    if timeStep % 10 == 0:                
        #dolfin.plot(dolfin.sqrt(dolfin.inner(NSDomain._solver.u1,NSDomain._solver.u1)),key='vNorm')         
        vFile << (eulerian._EulerianSolver__solver.u1, T)
        pFile << (eulerian._EulerianSolver__solver.p1, T)
        wFile << (eulerian._EulerianSolver__solver.vorticity(), T)

    print "Step\t\t\t: %g" % timeStep
    print "T\t\t\t: %g" % T
    print "Difference in vel. norm : %g" % diff_u    
    print "----------------------------------------\n" 
    
    if diff_u <= 1e-10:
        break
    
