"""

Test the navier-stokes modules: _navierStokes.py

"""

import pHyFlow

# External modules
import dolfin
import time
import numpy as np
import pylab as py
py.ion()


nTimeSteps = 50

# Geometry Position
cmGlobal = np.array([0.,0.])
thetaLocal = np.deg2rad(0.)

# Fluid Parameters
Re = 550.
U = 1.0
D = 2.0
nu   = (U*D)/(Re) # Unit Geometry with Re = 1000.

# CFL Number
cfl = 0.95
uMax = U*2.5

# Probe Grid Parameters
x0,y0,Lx,Ly,nx,ny = -1.0,-1.0,2.0,2.0,10,10    
#-----------------------------------------------------------------------------




#-----------------------------------------------------------------------------   
# Mmesh and boundary mesh files

# Export to files
meshFileName = './geometry/cylinder2D_hybrid_nVertices-13k_delQuad_nVertices-47k_pressureOutlet.xml.gz'
boundaryDomainsFileName = './geometry/cylinder2D_hybrid_nVertices-13k_delQuad_nVertices-47k_pressureOutlet_facet_region.xml.gz'

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
eulerian = pHyFlow.eulerian.EulerianSolver(geometry,probeGrid,uMax,nu,cfl)
                                     
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
vxBoundary, vyBoundary = 1.0*np.ones(xB.shape), np.zeros(yB.shape)
cmGlobalNew, thetaGlobalNew = cmGlobal,thetaLocal
cmDotGlobal, thetaGlobal = np.array([0.,0.]), 0.

saveDir = './data/'
vFile = dolfin.File(saveDir + "velocity.pvd", "compressed")
pFile = dolfin.File(saveDir + "pressure.pvd", "compressed")
wFile = dolfin.File(saveDir + "vorticity.pvd", "compressed")

solver = eulerian._EulerianSolver__solver

# Forces
CL = np.zeros(nTimeSteps+1)
CD = np.zeros(nTimeSteps+1)
CDfric = np.zeros(nTimeSteps+1)
CDpres = np.zeros(nTimeSteps+1)

for timeStep in range(1,nTimeSteps+1):

    startTime = time.time()
    
    # Current time
    T = timeStep*eulerian.deltaT

    # Evolve eulerian solution
    eulerian.evolve(vxBoundary,vyBoundary,cmGlobalNew,thetaGlobalNew,
                    cmDotGlobal,thetaGlobal)
    
    # Determine the norm on velocity
    diff_u = dolfin.norm(solver.u1) - dolfin.norm(solver.u0)            
    
    # Calculate the body forces
    CD[timeStep], CL[timeStep] = eulerian.Forces() / (0.5*U*U*D)
        
    # Pressure drag
    CDpres[timeStep] = eulerian.PressureForces()[0] / (0.5*U*U*D)
    # Frictional Drag
    CDfric[timeStep] = eulerian.FrictionalForces()[0] / (0.5*U*U*D)

    # Export results
    if timeStep % 1 == 0:                
        #dolfin.plot(dolfin.sqrt(dolfin.inner(NSDomain._solver.u1,NSDomain._solver.u1)),key='vNorm')         
        vFile << (eulerian._EulerianSolver__solver.u1, T)
        pFile << (eulerian._EulerianSolver__solver.p1, T)
        wFile << (eulerian._EulerianSolver__solver.vorticity(), T)

    # Print info
    print "Step\t\t\t: %g" % timeStep
    print "T\t\t\t: %g" % T
    print "Step duration\t\t: %g" % (time.time() - startTime)
    print "Difference in vel. norm : %g" % diff_u    
    print "----------------------------------------\n"
    
#------------------------------------------------------------------------------    
# Plot results
# Latex Text
py.rc('text', usetex=True)
py.rc('font', family='serif')


py.figure()
py.plot(np.arange(0,nTimeSteps+1)*eulerian.deltaT,CL)
py.xlabel(r'$t [-]$')
py.ylabel(r'$C_L$')
py.axis([0., (nTimeSteps+1)*eulerian.deltaT, -0.2,0.2])
py.grid()

py.figure()
py.plot(np.arange(0,nTimeSteps+1)*eulerian.deltaT,CD,'b-',label='Total')
py.plot(np.arange(0,nTimeSteps+1)*eulerian.deltaT,CDfric,'b-.',label='Friction')
py.plot(np.arange(0,nTimeSteps+1)*eulerian.deltaT,CDpres,'b--',label='Pressure')
py.xlabel(r'$t [-]$')
py.ylabel(r'$C_D$')
py.axis([0., (nTimeSteps+1)*eulerian.deltaT, -5,5])
py.legend()
py.grid()
    
