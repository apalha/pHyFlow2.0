"""

Hybrid test : ex01_lambOseen_test
=================================

Goal
----
 - test the coupling of navier-stokes and vortex blobs
     without any panel body.

 - test the evolution of lamb-oseen vortex core     
 
:First added:   2014-02-27
:Last updated:  2014-03-03
:Copyright:     Copyright (C) 2014 Lento Manickathan, **pHyFlow**
:Licence:       GNU GPL version 3 or any later version
"""

# External modules
import dolfin
import numpy as np
import pylab as py
import time
py.ion()


# pHyFlow
import pHyFlow

#-----------------------------------------------------------------------------    
# Global parameters

# Simulation parameters
nTimeSteps = 100

# Fluid Parameters
Re = 100. # Reynolds number
vInf = np.array([0.,0.]) # Free stream velocity


# Define the lamb-oseen
Gamma = 2.0*np.pi # Total circulation of the lamb-oseen
nu = Gamma / (2.0*np.pi*Re) # Unit Geometry with Re = 1000.
tau = 0.5 # viscous time
tStart = tau/nu
pVortex = np.array([0.,0.])


# Exact Function : exact vorticity field
def wLambOseen(x,y,t):
    
    return (Gamma / (4.0*np.pi*nu*t)) * np.exp( -  ((x-pVortex[0])*(x-pVortex[0]) +
                                                    (y-pVortex[1])*(y-pVortex[1]))/
                                                (4.0*nu*t))
                                                
# Exact Function : exact velocity field
def vLambOseen(x,y,t):
    # Radius
    r = np.sqrt((x-pVortex[0])*(x-pVortex[0]) + (y-pVortex[1])*(y-pVortex[1])) + np.spacing(1)
    # Circumferential velocity
    uTheta = (Gamma / (2.0*np.pi*r)) * (1.0 - np.exp( - (r*r)/(4.0*nu*t)))
    # Angle
    theta = np.arctan2((y-pVortex[1]),(x-pVortex[0]))
    # Return the cartesian velocity field
    return -np.sin(theta)*uTheta, np.cos(theta)*uTheta                                               
                                                
#-----------------------------------------------------------------------------    

#-----------------------------------------------------------------------------
# Setup the vortex blobs - Lagrangian domain

# Define blob parameters
overlap = 1.0               # overlap ratio
nBlobs = (256*256)    # number of blobs
h = 20.0/np.sqrt(nBlobs)    # blob spacing

# Estimate of delta T
#   with optimized c2 of 1/3.
#   delta T of convection should be a multiple of delta diffusion
deltaTc = float(np.around((1.0/5.0) * ( (1./3.)*h*h/nu),decimals=2))

# blob control parameters
blobControlParams = {'stepRedistribution':1,'stepPopulationControl':1,\
                     'gThresholdLocal':1e-8,'gThresholdGlobal':1e-8}


# Generate the vortex blobs uniformly distributed in the square [-10,10] x [-10,10]
x,y = np.meshgrid(np.linspace(-10.0,10.0,np.sqrt(nBlobs)),
                  np.linspace(-10.0,10.0,np.sqrt(nBlobs)))
x = x.flatten() + h*0.5
y = y.flatten() + h*0.5
g = h*h*wLambOseen(x,y,tStart)

wField = (x,y,g)

# generate the blobs
blobs = pHyFlow.blobs.Blobs(wField,vInf,nu,deltaTc,h,overlap,
                            blobControlParams=blobControlParams)

blobs.plotVelocity = True
# This problem does not have any panels
panels = None

# Setup the full vortex problem : blobs and others.
lagrangian = pHyFlow.lagrangian.LagrangianSolver(blobs,panels=None)

#------------------------------------------------------------------------------


#-----------------------------------------------------------------------------    
# Setup the navier-stokes problem - Eulerian domain

# Setup navier-stokes mesh and boundary Domains
#---------------------------

# Mesh parameters
N = 64# Number of mesh nodes (in x and y dir)

mesh = dolfin.UnitSquareMesh(N,N)
xy = mesh.coordinates()
xy = (xy-0.5)*2.0
mesh.coordinates()[:] = xy

# Sub domain for dirichlet boundary
class dirichletBoundary(dolfin.SubDomain):
    def inside(self,x,on_boundary):
        return on_boundary

# Define boundary domains    
boundaryDomains = dolfin.MeshFunction('size_t', mesh, mesh.topology().dim()-1)

# Mark the fluid
boundaryDomains.set_all(1)

# Mark no-slip
noslip = dirichletBoundary()
noslip.mark(boundaryDomains, pHyFlow.eulerian.eulerianOptions.ID_EXTERNAL_BOUNDARY)

# Export mesh and boundary mesh files

# Export to files
meshFilePath = './data/unitSquareMesh_%gx%g_mesh.xml.gz' % (N,N)
boundaryDomainsFilePath = './data/unitSquareMesh_%gx%g_boundaryDomains.xml.gz' % (N,N)
dolfin.File(meshFilePath) << mesh
dolfin.File(boundaryDomainsFilePath) << boundaryDomains
#---------------------------

# Define the geometry parameters
geometry = {'mesh' : meshFilePath,
            'boundaryDomains': boundaryDomainsFilePath,
            'cmGlobal': np.array([0.,0.]),
            'thetaLocal': 0.} # radians

# Probe Grid Parameters
probeL = np.array([3.0,3.0]) # guess
#origin = np.array([-1.0,-1.0])
probeN = np.int64(np.round(probeL / h))
probeL = probeN*h # to ensure the correct grid spacing
origin = -np.round(probeN*0.5)*h

probeGridParams = {'origin' : origin,
                   'L' : probeL,
                   'N' : probeN}          

# Solver parameters
cfl = 0.95 # CFL number
uMax = 1.5


# Initialize the navierStokes problem
eulerian = pHyFlow.eulerian.EulerianSolver(geometry,probeGridParams,
                                           uMax=uMax,nu=nu,cfl=cfl,
                                           deltaT=0.01)
  
# Concatenate all the navier-stokes domains
# TODO: temporary, will implement a multi-ns class
class MultiEulerianSolver(object):
    def __init__(self,subDomains):
        self.subDomains = subDomains
        
    def __getitem__(self,key):
        return self.subDomains[key]

    @property
    def nSubDomains(self):
        return len(self.subDomains)
        
    @property
    def subDomainKeys(self):
        return self.subDomains.keys()
        
    @property
    def t(self):
        return self.subDomains[self.subDomainKeys[0]].t

    @property
    def deltaT(self):
        return self.subDomains[self.subDomainKeys[0]].deltaT     

    @deltaT.setter
    def deltaT(self,deltaTNew):
        self.subDomains[self.subDomainKeys[0]].deltaT  = deltaTNew
     
    def getCoordinates(self):
        xyCoordinates = []
        for subDomain in self.subDomains.itervalues():
            xyCoordinates.append(subDomain.getCoordinates())
        
        return np.hstack(xyCoordinates)
        
    def getBoundaryCoordinates(self):
        xyBoundaryCoordinates = []
        for subDomain in self.subDomains.itervalues():
            xyBoundaryCoordinates.append(subDomain.getBoundaryCoordinates())
        
        return np.hstack(xyBoundaryCoordinates)   
        
    def getBoundaryVelocity(self):
        vxyBoundary = []
        for subDomain in self.subDomains.itervalues():
            vxyBoundary.append(subDomain.getBoundaryVelocity())
        
        return np.hstack(vxyBoundary)         
        
    def setVelocity(self,vx,vy):
        iS = 0
        for subDomain in self.subDomains.itervalues():
            iE = iS + subDomain.getCoordinates()[0].shape[0]
            subDomain.setVelocity(vx[iS:iE],vy[iS:iE])
            iS = np.copy(iE)

    def evolve(self,vx,vy,cmGlobalNew,thetaLocalNew,cmDotGlobalNew,thetaDotLocalNew):
        for subDomain in self.subDomains.itervalues():
            subDomain.evolve(vx,vy,cmGlobalNew,thetaLocalNew,cmDotGlobalNew,thetaDotLocalNew)
        
                                         
   
multiEulerian = MultiEulerianSolver(subDomains={'eulerian': eulerian})
        
                               
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Define the interpolation region

# Use stocks parameters
# ..[1] Stock, M., Gharakhani, A., & Stone, C. (2010). Modeling Rotor Wakes 
#       with a Hybrid OVERFLOW-Vortex Method on a GPU Cluster.
dBdry = 2*h + np.spacing(10e10)

interpolationRegions = {}

for subDomain in multiEulerian.subDomainKeys:

    # Get the limits of the navier-stokes boundary
    xMin = multiEulerian[subDomain].getBoundaryCoordinates()[0].min()
    xMax = multiEulerian[subDomain].getBoundaryCoordinates()[0].max()
    yMin = multiEulerian[subDomain].getBoundaryCoordinates()[1].min()
    yMax = multiEulerian[subDomain].getBoundaryCoordinates()[1].max()
    
    # Define the boundary polygon (closed loop)
    xBoundaryPolygon = np.array([xMin+dBdry, xMax-dBdry, xMax-dBdry,xMin+dBdry,xMin+dBdry])
    yBoundaryPolygon = np.array([yMin+dBdry, yMin+dBdry, yMax-dBdry,yMax-dBdry,yMin+dBdry])
    
    # Reshape the boundary polygon of point in polygon search
    xyBoundaryPolygon = np.vstack((xBoundaryPolygon,yBoundaryPolygon)).T
    
    interpolationRegions[subDomain] = {'surfacePolygon': None,
                                       'boundaryPolygon': xyBoundaryPolygon}


#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Setup the hybrid problem

# Define the coupling parameters
couplingParams={'adjustLagrangian':True,
                'adjustLagrangianAt':'start',
                'eulerianInitialConditions': 'lagrangian_field'}

interpolationParams={'algorithm':'scipy_griddata',
                     'method':'linear'}

# Initialize the hybrid problem
hybrid = pHyFlow.hybrid.HybridSolver(lagrangian=lagrangian, # either vortex-blobs or coupled vortex-panels
                                     multiEulerian=multiEulerian, # dictionary or class of multiple navier-stokes
                                     interpolationRegions=interpolationRegions, # the interpolation regions
                                     couplingParams=couplingParams,
                                     interpolationParams=interpolationParams) # coupling parameters
#-----------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Evolve

# Export the navier-stokes and blobs

# Define navier-stokes export files
nsVFileName = './data/Re%s/ns_velocity.pvd' % str(Re)
nsVFile = dolfin.File(nsVFileName)

nsWFileName = './data/Re%s/ns_vorticity.pvd' % str(Re)
nsWFile = dolfin.File(nsWFileName)

blobsFileName = './data/Re%s/blobs.pvd' % str(Re)
blobsFile = pHyFlow.IO.File(blobsFileName)


# Export initial data
nsVFile << (hybrid.multiEulerian['eulerian']._EulerianSolver__solver.u1, hybrid.multiEulerian.t)
nsWFile << (hybrid.multiEulerian['eulerian']._EulerianSolver__solver.vorticity(), hybrid.multiEulerian.t)
blobsFile << hybrid.lagrangian.blobs

#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Evolve

# new parameter
cmGlobalNew = np.array([0.,0.]) # same position as before
thetaLocalNew = 0. # same position as before
cmDotGlobalNew = np.array([0.,0.]) # zero velocity 
thetaDotLocalNew = 0. # zero-rotational velocity

# Post-processing parameters
numBlobsTime = np.zeros(nTimeSteps+1)
totalCirculationTime = np.zeros(nTimeSteps+1)
numBlobsTime[0] =  hybrid.lagrangian.blobs.numBlobs
totalCirculationTime[0] =  hybrid.lagrangian.blobs.g.sum()


for timeStep in range(1,nTimeSteps+1):
    
    startTime = time.time()
    # Evolve the coupled navier-stokes and vortex method
    hybrid.evolve(cmGlobalNew, thetaLocalNew, cmDotGlobalNew, thetaDotLocalNew)

    print "Time-step\t\t: %g" % timeStep
    print "T\t\t\t: %g" % hybrid.lagrangian.t
    print 'Time to evolve\t\t: %g'  % (time.time() - startTime)
    
    # Export initial data
    nsVFile << (hybrid.multiEulerian['eulerian']._EulerianSolver__solver.u1, hybrid.multiEulerian.t)
    nsWFile << (hybrid.multiEulerian['eulerian']._EulerianSolver__solver.vorticity(), hybrid.multiEulerian.t)
    blobsFile << hybrid.lagrangian.blobs
    
    numBlobsTime[timeStep] =  hybrid.lagrangian.blobs.numBlobs
    totalCirculationTime[timeStep] =  hybrid.lagrangian.blobs.g.sum()
    print 'Number of blobs\t\t: %d'  % numBlobsTime[timeStep]
    print 'Total circulation\t: %.8f'  % totalCirculationTime[timeStep]
    print '----------------------------------------\n'

#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
# Plot results

py.figure()
py.plot(np.arange(0,nTimeSteps+1)*hybrid.lagrangian.deltaT,numBlobsTime)
py.grid()
py.xlabel('t')

py.figure()
py.semilogy(np.arange(0,nTimeSteps+1)*hybrid.lagrangian.deltaT,np.abs(totalCirculationTime-totalCirculationTime[0]))
py.grid()
py.xlabel('t')

