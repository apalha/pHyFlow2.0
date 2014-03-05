"""

Hybrid test : ex02_cylinder_test
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
nTimeSteps = 10

# Fluid Parameters

# Double-monopole parameters
Re = 550. # Reynolds number
vInf = np.array([1.,0.]) # Free stream velocity
L = 2.0

# Compute fluid kinematic viscosity
nu = float(vInf[0]*L/Re)

# Define cylinder
Rext = 1.5
Rsurf = 1.0

#-----------------------------------------------------------------------------    

#-----------------------------------------------------------------------------
# Setup the vortex blobs - Lagrangian domain

# Define blob parameters
overlap = 1.0               # overlap ratio
nBlobs = (60*(2*Rext))**2# number of blobs
hBlob = (2*Rext)/np.sqrt(nBlobs)    # blob spacing

# Estimate of delta T
#   with optimized c2 of 1/3.
#   delta T of convection should be a multiple of delta diffusion
deltaTc = 0.01#float(np.round((1.0/nConvectionPerDiffusion) * ( (1./3.)*h*h/nu),decimals=3))

# blob control parameters
blobControlParams = {'stepRedistribution':1,'stepPopulationControl':1,\
                     'gThresholdLocal':1e-8,'gThresholdGlobal':1e-8}


# Generate the vortex blobs uniformly distributed in the square [-10,10] x [-10,10]
x,y = np.meshgrid(np.linspace(-Rext,Rext,np.sqrt(nBlobs)+1),
                  np.linspace(-Rext,Rext,np.sqrt(nBlobs)+1))
x = x.flatten()
y = y.flatten()
g = np.zeros(x.shape)

if np.sum(x==0) == 0 or np.sum(y==0) == 0:
    raise ValueError('vortex mesh should coincide with (0.,0.) for remeshing grid')

wField = (x,y,g)

# generate the blobs
blobs = pHyFlow.vortex.VortexBlobs(wField,vInf,nu,deltaTc,hBlob,overlap,
                                   blobControlParams=blobControlParams)

blobs.plotVelocity = True
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Setup the panels - Lagrangian domain
# This problem does not have any panels
nPanels = 100
dPanel = np.spacing(100)
theta   = np.linspace(np.pi,-np.pi,nPanels+1) # Panel polar angles, to define them.
dtheta  = theta[1]-theta[0] # Angle spacing
r       = (Rsurf + dPanel) / np.cos(dtheta/2.0) # Radial location of the panel end points

# Make the cylinder and append the parameters to a dictionary.
cylinderData = {'xPanel' : r*np.cos(theta - dtheta/2),
                'yPanel' : r*np.sin(theta - dtheta/2),
                'cmGlobal'   : np.array([0.,0.]),
                'thetaLocal' : 0.,
                'dPanel' : np.spacing(100)}
                
geometries = {'cylinder':cylinderData}                

panels = pHyFlow.panel.Panels(geometries=geometries)

#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Setup the full vortex problems

# Setup the full vortex problem : blobs and others.
vortex = pHyFlow.vortexPanel.VortexPanel(blobs,panels=panels,
                                         couplingParams={'panelStrengthUpdate':'constant'})

#------------------------------------------------------------------------------


#-----------------------------------------------------------------------------    
# Setup the navier-stokes problem - Eulerian domain

# Define the mesh

meshFilePath = './geometry/cylinder2D_nVertices-13k_delQuad_noBL_finerBody.xml.gz'
boundaryDomainsFilePath = './geometry/cylinder2D_nVertices-13k_delQuad_noBL_finerBody_facet_region.xml.gz'


# Define the geometry parameters
geometry = {'mesh' : meshFilePath,
            'boundaryDomains': boundaryDomainsFilePath,
            'cmGlobal': np.array([0.,0.]),
            'thetaLocal': 0.} # radians

# Probe Grid Parameters
probeL = np.array([3.0,3.0]) # guess
probeN = np.int64(np.round(probeL / hBlob))
probeL = probeN*hBlob # to ensure the correct grid spacing
origin = -np.round(probeN*0.5)*hBlob

# Collect all the probe parameters
probeGridParams = {'origin' : origin,
                   'L' : probeL,
                   'N' : probeN}          

# Solver parameters
cfl = 0.95 # CFL number
uMax = 2.5*vInf[0] # maximum fluid velocity

# Initialize the navierStokes problem
NSDomain = pHyFlow.navierStokes.NavierStokes(geometry,probeGridParams,
                                             uMax=uMax,nu=nu,cfl=cfl,
                                             deltaT=0.001)
  
# Concatenate all the navier-stokes domains
# TODO: temporary, will implement a multi-ns class
class MultiNavierStokes(object):
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
        
                                         
# Define the multiple naviers-stokes class   
multiNS = MultiNavierStokes(subDomains={'NSDomain': NSDomain})
        
                               
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Define the interpolation region

# Use stocks parameters
# ..[1] Stock, M., Gharakhani, A., & Stone, C. (2010). Modeling Rotor Wakes 
#       with a Hybrid OVERFLOW-Vortex Method on a GPU Cluster.
dBdry = Rsurf*0.1#2*hBlob  + np.spacing(10e10)
dSurf = (1.0 + 2.0) * hBlob

# Generate polar points (N = 100)
beta = np.linspace(-np.pi,np.pi,100)

# Define the boundary (exterior) polygon
xyBoundaryPolygon = np.vstack(( (Rext - dBdry)*np.cos(beta) + NSDomain.cmGlobal[0],
                                (Rext - dBdry)*np.sin(beta) + NSDomain.cmGlobal[1])).T
                                       
# Define the surface polygon
xySurfacePolygon =  np.vstack(( (Rsurf + dSurf)*np.cos(beta) + NSDomain.cmGlobal[0],
                                (Rsurf + dSurf)*np.sin(beta) + NSDomain.cmGlobal[1])).T
                               
interpolationRegions = {multiNS.subDomainKeys[0] : {'surfacePolygon': xySurfacePolygon,
                                                    'boundaryPolygon': xyBoundaryPolygon}}                      
                       

#-----------------------------------------------------------------------------
#
#-----------------------------------------------------------------------------
# Setup the hybrid problem

# Define the coupling parameters
couplingParams={'adjustLagrangian':True,
                'adjustLagrangianAt':'start',
                'eulerianInitialConditions': 'lagrangian_field'}

interpolationParams={'algorithm':'scipy_griddata',#'structured_probes',
                     'method':'linear'}

# Initialize the hybrid problem
hybrid = pHyFlow.hybrid.Hybrid(vortex=vortex, # either vortex-blobs or coupled vortex-panels
                               navierStokes=multiNS, # dictionary or class of multiple navier-stokes
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
nsVFile << (hybrid.navierStokes['NSDomain']._NavierStokes__solver.u1, hybrid.navierStokes.t)
nsWFile << (hybrid.navierStokes['NSDomain']._NavierStokes__solver.vorticity(), hybrid.navierStokes.t)
#blobsFile << blobs

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
numBlobsTime[0] =  hybrid.vortex.blobs.numBlobs
totalCirculationTime[0] =  hybrid.vortex.blobs.g.sum()


for timeStep in range(1,nTimeSteps+1):
    
    
    
    startTime = time.time()
    # Evolve the coupled navier-stokes and vortex method
    hybrid.evolve(cmGlobalNew, thetaLocalNew, cmDotGlobalNew, thetaDotLocalNew)

    print "\nThe current time\t\t: %g" % (timeStep*hybrid.vortex.deltaT)
    print "Time-step\t\t\t: %g" % timeStep
    print 'Time to evolve\t\t\t: %g'  % (time.time() - startTime)
    
    # Export initial data
    if timeStep < 100:
            nsVFile << (hybrid.navierStokes['NSDomain']._NavierStokes__solver.u1, hybrid.navierStokes.t)
            nsWFile << (hybrid.navierStokes['NSDomain']._NavierStokes__solver.vorticity(), hybrid.navierStokes.t)
            blobsFile << hybrid.vortex.blobs
    else:
        if timeStep % 50 == 0:
            nsVFile << (hybrid.navierStokes['NSDomain']._NavierStokes__solver.u1, hybrid.navierStokes.t)
            nsWFile << (hybrid.navierStokes['NSDomain']._NavierStokes__solver.vorticity(), hybrid.navierStokes.t)
            blobsFile << hybrid.vortex.blobs
    
    numBlobsTime[timeStep] =  hybrid.vortex.blobs.numBlobs
    totalCirculationTime[timeStep] =  hybrid.vortex.blobs.g.sum()
    print 'Number of blobs\t\t\t: %d'  % numBlobsTime[timeStep]
    print 'Total circulation\t\t: %.8f'  % totalCirculationTime[timeStep]
    print '--------------------------------------------\n'

#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
# Plot results
# Latex Text
py.rc('text', usetex=True)
py.rc('font', family='serif')


py.figure()
py.plot(np.arange(0,nTimeSteps+1)*blobs.deltaTc,numBlobsTime)
py.xlabel(r'$t [-]$')
py.ylabel(r'$N_{blobs}$')
py.grid()

py.figure()
py.semilogy(np.arange(0,nTimeSteps+1)*hybrid.vortex.deltaT,np.abs(totalCirculationTime-totalCirculationTime[0]))
py.xlabel(r'$t [-]$')
py.ylabel(r'$\Delta \Gamma$')
py.grid()
