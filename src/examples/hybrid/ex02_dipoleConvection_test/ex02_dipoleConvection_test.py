"""

Hybrid test : ex02_dipoleConvection_test
=================================

Goal
----
 - test the coupling of navier-stokes and vortex blobs
     without any panel body.

 - test the evolution of lamb-oseen vortex core     

:First added:   2014-03-03
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
nTimeSteps = 50

# Fluid Parameters

# Double-monopole parameters
Re = 625. # Reynolds number
vInf = np.array([0.,0.]) # Free stream velocity
U = 1.0
L = 1.0
Tend = 1.0

# Compute fluid kinematic viscosity
nu = U*L/Re         

# Initial conditions for a double-Monopole
# Parameters
x1,y1   = 0.1, 1. # Initial core location. `Positive Core`
x2,y2   = -0.1,1. # Initial core location. `Negative Core`
r0      = 0.1      # Half core radius
omega_e = 299.528385375226# Maximum single core voriticty

# Exact Function
def wExactFunction(x,y):
    rr1 = (x-x1)**2 + (y-y1)**2 # Radius 1
    rr2 = (x-x2)**2 + (y-y2)**2 # Radius 2
    return omega_e*(1- rr1/(r0**2))*np.exp(-rr1/(r0**2)) \
          -omega_e*(1- rr2/(r0**2))*np.exp(-rr2/(r0**2))
                                  
#-----------------------------------------------------------------------------    

#-----------------------------------------------------------------------------
# Setup the vortex blobs - Lagrangian domain

# Define blob parameters
overlap = 1.0               # overlap ratio
nBlobs = 300*300# number of blobs
h = 2.0/np.sqrt(nBlobs)    # blob spacing

# Estimate of delta T
#   with optimized c2 of 1/3.
#   delta T of convection should be a multiple of delta diffusion
#nConvectionPerDiffusion = 1.
deltaTc = 0.002#float(np.round((1.0/nConvectionPerDiffusion) * ( (1./3.)*h*h/nu),decimals=3))

# blob control parameters
blobControlParams = {'stepRedistribution':1,'stepPopulationControl':1,\
                     'gThresholdLocal':1e-8,'gThresholdGlobal':1e-8}


# Generate the vortex blobs uniformly distributed in the square [-10,10] x [-10,10]
x,y = np.meshgrid(np.linspace(-1.0,1.0,np.sqrt(nBlobs)+1),
                  np.linspace(y1-1.0,y1+1.0,np.sqrt(nBlobs)+1))
x = x.flatten()
y = y.flatten()
g = wExactFunction(x,y) * (h*h)

if np.sum(x==0) == 0 or np.sum(y==0) == 0:
    raise ValueError('vortex mesh should coincide with (0.,0.) for remeshing grid')

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

## Mesh parameters
#N = 100# Number of mesh nodes (in x and y dir)
#
#mesh = dolfin.UnitSquareMesh(N,N)
#xy = mesh.coordinates()
#xy = (xy - 0.5)*0.5
#mesh.coordinates()[:] = xy
#
## Sub domain for dirichlet boundary
#class dirichletBoundary(dolfin.SubDomain):
#    def inside(self,x,on_boundary):
#        return on_boundary
#
## Define boundary domains    
#boundaryDomains = dolfin.MeshFunction('size_t', mesh, mesh.topology().dim()-1)
#
## Mark the fluid
#boundaryDomains.set_all(1)
#
## Mark no-slip
#noslip = dirichletBoundary()
#noslip.mark(boundaryDomains, pHyFlow.navierStokes.nsOptions.ID_EXTERNAL_BOUNDARY)

# Export mesh and boundary mesh files

# Export to files
#meshFilePath = './data/unitSquareMesh_%gx%g_mesh.xml.gz' % (N,N)
#boundaryDomainsFilePath = './data/unitSquareMesh_%gx%g_boundaryDomains.xml.gz' % (N,N)
#dolfin.File(meshFilePath) << mesh
#dolfin.File(boundaryDomainsFilePath) << boundaryDomains

meshFilePath = './geometry/dipoleConvection_Re-625_vertices-2k_MPI.xml.gz'
boundaryDomainsFilePath = './geometry/dipoleConvection_Re-625_vertices-2k_MPI_facet_region.xml.gz'
#---------------------------


# Define the geometry parameters
geometry = {'mesh' : meshFilePath,
            'boundaryDomains': boundaryDomainsFilePath,
            'cmGlobal': np.array([0.,0.]),
            'thetaLocal': 0.} # radians

# Probe Grid Parameters
probeL = np.array([1.0,0.5]) # guess
origin = np.array([-0.5,-0.25])
probeN = np.int64(np.round(probeL / h))
probeL = probeN*h # to ensure the correct grid spacing
origin = -np.round(probeN*0.5)*h

probeGridParams = {'origin' : origin,
                   'L' : probeL,
                   'N' : probeN}          

# Solver parameters
cfl = 0.95 # CFL number
uMax = 12.0


# Initialize the navierStokes problem
eulerian = pHyFlow.eulerian.EulerianSolver(geometry,probeGridParams,
                                           uMax=uMax,nu=nu,cfl=cfl,
                                           deltaT=5e-5)
  
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
        
                                         
# Initialize the eulerian solver   
multiEulerian = MultiEulerianSolver(subDomains={'eulerian': eulerian})
        
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Define the interpolation region

# Use stocks parameters
# ..[1] Stock, M., Gharakhani, A., & Stone, C. (2010). Modeling Rotor Wakes 
#       with a Hybrid OVERFLOW-Vortex Method on a GPU Cluster.
dBdry = 2*h  + np.spacing(10e10) # line should not coincide with remeshing grid

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
hybrid = pHyFlow.hybrid.HybridSolver(lagrangian=lagrangian,
                                     multiEulerian=multiEulerian, # dictionary or class of multiple navier-stokes
                                     interpolationRegions=interpolationRegions, # the interpolation regions
                                     couplingParams=couplingParams,
                                     interpolationParams=interpolationParams) # coupling parameters
#-----------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Evolve

# Export the navier-stokes and blobs

# Define navier-stokes export files
nsVFileName = './data/ns_velocity.pvd'
nsVFile = dolfin.File(nsVFileName)

nsWFileName = './data/ns_vorticity.pvd'
nsWFile = dolfin.File(nsWFileName)

blobsFileName = './data/blobs.pvd'
blobsFile = pHyFlow.IO.File(blobsFileName)


# Export initial data
nsVFile << (hybrid.multiEulerian['eulerian']._EulerianSolver__solver.u1, hybrid.multiEulerian.t)
blobsFile << blobs

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
    print "T\t\t\t: %g" % (timeStep*hybrid.lagrangian.deltaT)
    print 'Time to evolve\t\t: %g'  % (time.time() - startTime)
    
    # Export initial data
    if timeStep % 1 == 0:
        nsVFile << (hybrid.multiEulerian['eulerian']._EulerianSolver__solver.u1, hybrid.multiEulerian.t)
        nsWFile << (hybrid.multiEulerian['eulerian']._EulerianSolver__solver.vorticity(), hybrid.multiEulerian.t)
        blobsFile << hybrid.lagrangian.blobs
    
    numBlobsTime[timeStep] =  hybrid.lagrangian.blobs.numBlobs
    totalCirculationTime[timeStep] =  hybrid.lagrangian.blobs.g.sum()
    print 'Number of blobs\t\t: %g'  % numBlobsTime[timeStep]
    print 'Total circulation\t: %g'  % totalCirculationTime[timeStep]
    print '--------------------------------------------\n'

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

