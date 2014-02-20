"""

Test the navier-stokes modules: _navierStokes.py

"""

import pHyFlow

# External modules
import time
import dolfin
import numpy as np
import pylab as py
py.ion()

#def initialize_from_parameters():
#    """
#    Test using Lid-Driven cavity Flow
#    """
#    pass
 

#-----------------------------------------------------------------------------    

# Mesh Data
N = 40 
#mesh, boundaryDomains = lidDrivenCavityFlow(N)

# Geometry Position
cmGlobal = np.array([0.,0.])
thetaGlobal = np.deg2rad(0.)

# Fluid Parameters
Re = 1000.
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
noslip.mark(boundaryDomains, pHyFlow.navierStokes.noSlipID) # No-slip ID: 2

# Mark exterior I.D
lid = Lid()
lid.mark(boundaryDomains, pHyFlow.navierStokes.extID) # extID : 3  
  
#-----------------------------------------------------------------------------   


#-----------------------------------------------------------------------------   
# Export mesh and boundary mesh files

# Export to files
meshFileName = 'LDCF_%gx%g_mesh.xml.gz' % (N,N)
boundaryMeshFileName = 'LDCF_%gx%g_boundaryDomains.xml.gz' % (N,N)
dolfin.File(meshFileName) << mesh
dolfin.File(boundaryMeshFileName) << boundaryDomains
#-----------------------------------------------------------------------------   


#------------------------------------------------------------------------------
# Lid-Driven cavity flow test case

# Initialize Navier-Stokes problem
NSDomain = pHyFlow.navierStokes.NavierStokes(mesh=meshFileName,
                                     boundaryDomains=boundaryMeshFileName,
                                     cmGlobal=cmGlobal,thetaGlobal=thetaGlobal,
                                     uMax=uMax,nu=nu,cfl=cfl,origin=np.array([x0,y0]),
                                     L=np.array([Lx,Ly]),N=np.array([nx,ny]))
                                     
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
# Plot Velocity vectors - ONLY FOR TESTING (else need to output to pvd)
vx,vy = np.random.rand(2,x.shape[0])
NSDomain.setVelocity(vx,vy)
dolfin.plot(NSDomain._NavierStokes__solver.u1, title='ONLY FOR TESTING - Random velocity field')
dolfin.interactive()

#------------------------------------------------------------------------------
# Set velocity
vx,vy = np.zeros((2,x.shape[0]))
NSDomain.setVelocity(vx,vy)
# PLlot ONLY FOR TESTING (else need to output to pvd)
dolfin.plot(NSDomain._NavierStokes__solver.u1, title='ONLY FOR TESTING - Changed to zero velocity field')
dolfin.interactive()
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
vxBoundary, vyBoundary = -1.0*np.ones(xB.shape), np.zeros(yB.shape)
cmGlobalNew, thetaGlobalNew = None,None
cmDotGlobal, thetaGlobal = None, None

saveDir = './data/'
vFile = dolfin.File(saveDir + "velocity.pvd", "compressed")
pFile = dolfin.File(saveDir + "pressure.pvd", "compressed")
wFile = dolfin.File(saveDir + "vorticity.pvd", "compressed")

solver = NSDomain._NavierStokes__solver
for i in range(50000):

    T = i*NSDomain.dtMax

    NSDomain.evolve(vxBoundary,vyBoundary,cmGlobalNew,thetaGlobalNew,
                    cmDotGlobal,thetaGlobal)
    
    diff_u = dolfin.norm(solver.u1) - dolfin.norm(solver.u0)            

    if i % 200 == 0:                
        #dolfin.plot(dolfin.sqrt(dolfin.inner(NSDomain._solver.u1,NSDomain._solver.u1)),key='vNorm')         
        vFile << (NSDomain._NavierStokes__solver.u1, T)
        pFile << (NSDomain._NavierStokes__solver.p1, T)
        wFile << (NSDomain._NavierStokes__solver.vorticity(), T)

        print "Step = %g\tT = %g" % (i, T)
        print "Difference in velocity : %g" % diff_u    
    
    if diff_u <= 1e-10:
        break
    