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

for timeStep in xrange(1,nTimeSteps+1):

    startTime = time.time()
    
    # Current time
    T = timeStep*eulerian.deltaT

    eulerian.evolve(vxBoundary,vyBoundary,cmGlobalNew,thetaGlobalNew,
                    cmDotGlobal,thetaGlobal)
    
    # Determine the norm on velocity
    diff_u = dolfin.norm(solver.u1) - dolfin.norm(solver.u0)            
    
    # Calculate the body forces
    CD[timeStep], CL[timeStep] = eulerian.Forces() / (0.5*uInf*uInf*L)
        
    # Pressure drag
    CDpres[timeStep] = eulerian.PressureForces()[0] / (0.5*uInf*uInf*L)
    # Frictional Drag
    CDfric[timeStep] = eulerian.FrictionalForces()[0] / (0.5*uInf*uInf*L)           

    if timeStep % 5 == 0:                
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
py.axis([0., (nTimeSteps+1)*eulerian.deltaT, -5,5])
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
       