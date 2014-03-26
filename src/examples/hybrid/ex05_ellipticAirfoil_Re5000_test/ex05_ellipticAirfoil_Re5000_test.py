"""

Hybrid test : ex05_airfoil_NACA0012_test
=========================================

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
Re = 5000. # Reynolds number
vInf = np.array([1.,0.]) # Free stream velocity

# Airfoil parameters
L = 1.0
t = 0.1
cmGlobal = np.array([0.,0.]) # location
cmLocal =  np.array([0.25,0.]) # location
thetaLocal = np.deg2rad(-20.)

xMeshBounds = np.array([-0.5,1.0])
yMeshBounds = np.array([-0.5,0.5])

# Compute fluid kinematic viscosity
nu = float(vInf[0]*L/Re)
#-----------------------------------------------------------------------------    

##-----------------------------------------------------------------------------
## Setup the vortex blobs - Lagrangian domain
#
# Define blob parameters
overlap = 1.0               # overlap ratio
nBlobs = (400*1.5*1)**2    # number of blobs
hBlob = (1.5)/np.sqrt(nBlobs)    # blob spacing
print "hBlob: %g" % hBlob
# Estimate of delta T
#   with optimized c2 of 1/3.
#   delta T of convection should be a multiple of delta diffusion
#nConvectionPerDiffusion = 1
deltaTc = 0.0002#float(np.round((1.0/nConvectionPerDiffusion) * ( (1./3.)*hBlob*hBlob/nu),decimals=4))

# blob control parameters
blobControlParams = {'stepRedistribution':1,'stepPopulationControl':1,\
                     'gThresholdLocal':1e-8,'gThresholdGlobal':1e-8}
#diffusionParams={'method':'regrid_wee','c2':1./6.}

# Generate the vortex blobs uniformly distributed in the square [-10,10] x [-10,10]
xi = np.arange(0,np.int64(400*1.5))*hBlob - np.floor((400*1.5)/2.0)*hBlob
yi = np.arange(0,np.int64(400*1.0))*hBlob - np.floor((400*1.0)/2.0)*hBlob

x,y = np.meshgrid(xi,yi)
x = x.flatten()
y = y.flatten()
g = np.zeros(x.shape)

if np.sum(x==0) == 0 or np.sum(y==0) == 0:
    raise ValueError('vortex mesh should coincide with (0.,0.) for remeshing grid')

wField = (x,y,g)

nuBlob = 0.

# generate the blobs
blobs = pHyFlow.blobs.Blobs(wField,vInf,nuBlob,deltaTc,hBlob,overlap,
                            blobControlParams=blobControlParams)

print "stepDiffusion : %g" % blobs.stepDiffusion
print "deltaTc : %g" % blobs.deltaTc
blobs.plotVelocity = True
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Setup the panels - Lagrangian domain
# This problem does not have any panels
nPanels = 100
dPanel = np.spacing(100)

# Generate sin airfoil mapping
xi = np.linspace(0,1,np.ceil(0.5*nPanels+1))
xxi = 0.5*L*(1.0 - np.cos(np.pi*xi))

#theta = np.linspace(np.pi,-np.pi,nPanels+1) # Panel polar angles
# Sin mapping of thetas
#theta = np.pi - 2*np.pi*xxi 


# Determine panel parameters
theta = np.linspace(np.pi,-np.pi,nPanels+1) # Panel polar angles
dtheta = theta[1]-theta[0]

# Ellipse parameters
a = L/2.0
b = t/2.0

# Panel Coordinates in cartesian coordinates
xPanel = a*np.cos(theta - dtheta/2) + cmLocal[0]
yPanel = b*np.sin(theta - dtheta/2) + cmLocal[1]   

# Make the cylinder and append the parameters to a dictionary.
panelData = {'xPanel' : xPanel,
            'yPanel' : yPanel,
            'cmGlobal'   : cmGlobal,
            'thetaLocal' : thetaLocal,
            'dPanel' : dPanel}
            
geometries = {'ellipseAirfoil':panelData}                

problemType='prescribed'

# Initialize the panels
panels = pHyFlow.panels.Panels(geometries=geometries,
                               panelKernel='csv',
                               problemType=problemType)

#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Setup the full vortex problems

# Setup the full vortex problem : blobs and others.
lagrangian = pHyFlow.lagrangian.LagrangianSolver(blobs,panels=panels,
                                                 couplingParams={'panelStrengthUpdate':'constant'})
#------------------------------------------------------------------------------


#-----------------------------------------------------------------------------
# Define the interpolation region

# Use stocks parameters
# ..[1] Stock, M., Gharakhani, A., & Stone, C. (2010). Modeling Rotor Wakes 
#       with a Hybrid OVERFLOW-Vortex Method on a GPU Cluster.
#dBdry = 0.4*0.1#2*hBlob  + np.spacing(10e10)
dSurf = 3.0*hBlob
#dBdry = 2.0*hBlob

# This problem does not have any panels
nSurfacePoints = 100

# Important : counter-clockwise
theta = np.linspace(-np.pi,np.pi,nSurfacePoints+1) # Panel polar angles

# Ellipse parameters
a = L/2.0 + dSurf
b = t/2.0 + dSurf

# Panel Coordinates in cartesian coordinates
x = a*np.cos(theta) + cmLocal[0]
y = b*np.sin(theta) + cmLocal[1]   

xySurfacePolygon = np.vstack((x,y)).T
                                       
# Define the surface polygon
nPoints = 100
# Determine panel parameters
theta = np.linspace(-np.pi,np.pi,nPoints+1) # Panel polar angles

# Ellipse parameters
a = 0.625*0.9
b = 0.375*0.9
#a = 0.625 - 2.0*hBlob
#b = 0.375 - 2.0*hBlob

# Panel Coordinates in cartesian coordinates
x = a*np.cos(theta) + cmLocal[0]
y = b*np.sin(theta) + cmLocal[1]                                       

# Boundary polygon coordinates
xyBoundaryPolygon =  np.vstack((x,y)).T

                                                          
#-----------------------------------------------------------------------------


#-----------------------------------------------------------------------------    
# Setup the navier-stokes problem - Eulerian domain
#
# Define the mesh

meshFilePath = './geometry/ellipticAirfoil_nVertices-59k.xml.gz'
boundaryDomainsFilePath = './geometry/ellipticAirfoil_nVertices-59k_facet_region.xml.gz'

xMeshBounds = np.array([-0.5,1.5])
yMeshBounds = np.array
# Define the geometry parameters
geometry = {'mesh' : meshFilePath,
            'boundaryDomains': boundaryDomainsFilePath,
            'cmGlobal': cmGlobal,
            'thetaLocal': thetaLocal} # radians

# Probe Grid Parameters
probeL = np.array([1.5,1.0]) # guess
probeN = np.int64(np.round(probeL / hBlob)) + 1 
probeL = (probeN-1)*hBlob # to ensure the correct grid spacing
origin = -np.round(probeN*0.5)*hBlob + cmLocal


# Boundary polygon coordinates
xyCirculationProbes = xyBoundaryPolygon.T.copy()

# Collect all the probe parameters
probeGridParams = {'origin' : origin,
                   'L' : probeL,
                   'N' : probeN,
                   'circulationProbes':xyCirculationProbes} # Important : counter-clockwise

# Solver parameters
cfl = 1.0 # CFL number
uMax = 3.0 # maximum fluid velocity

# Initialize the navierStokes problem
eulerian = pHyFlow.eulerian.EulerianSolver(geometry,probeGridParams,
                                           uMax=uMax,nu=nu,cfl=cfl,
                                           deltaT=0.0001)

eulerian.plotPressure = True                                           
eulerian.plotVelocity = True

#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# multi-eulerian class

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
            
    def gTotalInside(self):
        gTotalInsideList = []
        for subDomain in self.subDomains.itervalues():
            gTotalInsideList.append(subDomain.gTotalInside())
            
        return np.array(gTotalInsideList)
        
                                         
# Define the multiple naviers-stokes class   
multiEulerian = MultiEulerianSolver(subDomains={'eulerian': eulerian})

# Only one cylinder                               
interpolationRegions = {multiEulerian.subDomainKeys[0] : {'surfacePolygon': xySurfacePolygon,
                                                          'boundaryPolygon': xyBoundaryPolygon}}        
 
#-----------------------------------------------------------------------------


#-----------------------------------------------------------------------------
# Setup the hybrid problem

# Define the coupling parameters
couplingParams={'adjustLagrangian':True,
                'adjustLagrangianAt':'end',
                'eulerianInitialConditions': 'eulerian_field',
                'conserveCirculation':True}#'lagrangian_field'}

interpolationParams={'algorithm':'structuredProbes_manual',#structuredProbes_manual',#'structured_probes',
                     'method':'linear'}

# Initialize the hybrid problem
hybrid = pHyFlow.hybrid.HybridSolver(lagrangian=lagrangian, # either vortex-blobs or coupled vortex-panels
                                     multiEulerian=multiEulerian, # dictionary or class of multiple navier-stokes
                                     interpolationRegions=interpolationRegions, # the interpolation regions
                                     couplingParams=couplingParams,
                                     interpolationParams=interpolationParams) # coupling parameters
                                     
## Ensure zero velocity at no-slip boundary
#solver = eulerian._EulerianSolver__solver
#iX,iY = pHyFlow.eulerian.base.boundary.vectorDOF_boundaryIndex(solver.mesh,
#                                                               solver.boundaryDomains,
#                                                               pHyFlow.eulerian.eulerianOptions.ID_NOSLIP_BOUNDARY,
#                                                               2)
## Set vx and vy at no-slip boundary to zero
#solver.u1.vector()[iX] = 0.
#solver.u1.vector()[iY] = 0.

print "nEulerianSubSteps: %g " % hybrid.nEulerianSubSteps
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Determine lagrangian induced velocity

meshL = dolfin.RectangleMesh(-0.5,-0.5,1.0,0.5,100,100)

VL = dolfin.VectorFunctionSpace(meshL,'CG',1)
xL,yL = np.copy(VL.dofmap().tabulate_all_coordinates(meshL).reshape(-1,2).T[:,::2])

uLagrangian = dolfin.Function(VL)

#-----------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Export results

# Export the navier-stokes and blobs

# Post-proc save dir
postProcSaveDir = './data/'
postProcStep = 50

# Define eulerian export files
eulerianFile = pHyFlow.IO.File(postProcSaveDir+'eulerian.pvd')

# Define the blob export file
blobsFile = pHyFlow.IO.File(postProcSaveDir+'blobs.pvd')

# Define the panel export file
panelFile = pHyFlow.IO.File(postProcSaveDir+'panels.pvd')

lagrangianFieldFile = dolfin.File(postProcSaveDir+'lagrangianField.pvd')

# Export initial data
eulerianFile << hybrid.multiEulerian['eulerian']
panelFile << hybrid.lagrangian.panels
blobsFile << hybrid.lagrangian.blobs

vxL, vyL = hybrid.lagrangian.evaluateVelocity(xL,yL)
uLagrangian.vector()[:] = np.vstack((vxL, vyL)).T.flatten()

lagrangianFieldFile << (uLagrangian,lagrangian.t)



#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Evolve

# new parameter
cmGlobalNew = np.array([0.,0.]) # same position as before
thetaLocalNew = 0. # same position as before
cmDotGlobalNew = np.array([0.,0.]) # zero velocity 
thetaDotLocalNew = 0. # zero-rotational velocity

# Post-processing parameters
numBlobsTime = np.zeros((nTimeSteps/postProcStep)+1)
totalCirculationTime = np.zeros((nTimeSteps/postProcStep)+1)
totalPanelCirculationTime = np.zeros((nTimeSteps/postProcStep)+1)
totalBlobCirculationTime = np.zeros((nTimeSteps/postProcStep)+1)

numBlobsTime[0] =  hybrid.lagrangian.blobs.numBlobs
totalCirculationTime[0] =  hybrid.lagrangian.gTotal
totalPanelCirculationTime[0] =  hybrid.lagrangian.panels.gTotal[0]
totalBlobCirculationTime[0] = hybrid.lagrangian.blobs.g.sum()
print "totalCirculation: %g" % totalCirculationTime[0]
# Forces
CL = np.zeros((nTimeSteps/postProcStep)+1)
CD = np.zeros((nTimeSteps/postProcStep)+1)
CDfric = np.zeros((nTimeSteps/postProcStep)+1)
CDpres = np.zeros((nTimeSteps/postProcStep)+1)

# Calculate the body forces
CD[0], CL[0] = eulerian.Forces() / (0.5*vInf[0]*vInf[0]*L)
# Pressure drag
CDpres[0] = eulerian.PressureForces()[0] / (0.5*vInf[0]*vInf[0]*L)
# Frictional Drag
CDfric[0] = eulerian.FrictionalForces()[0] / (0.5*vInf[0]*vInf[0]*L)

for timeStep in range(1,nTimeSteps+1):
        
    # timer
    startTime = time.time()
    
    # Current time
    T = timeStep*hybrid.lagrangian.deltaT
    
    # Evolve the coupled navier-stokes and vortex method
    hybrid.evolve(cmGlobalNew, thetaLocalNew, cmDotGlobalNew, thetaDotLocalNew)
    
    # Print info
    endTime = (time.time() - startTime)
    timeLeft = endTime*(nTimeSteps+1-timeStep) 
    print "\nTime-step\t\t\t: %g" % timeStep
    print "Current time\t\t\t: %g" % (timeStep*hybrid.lagrangian.deltaT)        
    print 'Time to evolve\t\t\t: %g'  % (endTime)
    print 'End time\t\t\t: %g sec, %g hrs, %g days.' % (timeLeft, timeLeft/(60.*60.), timeLeft/(60.*60.*24.))
    print 'Current local time\t\t: %s' % time.strftime("%Y-%m-%d, %H:%M:%S", time.localtime())
    # Export/calculate data every postProcStep

    if timeStep < 5:
        eulerianFile << hybrid.multiEulerian['eulerian']
        panelFile << hybrid.lagrangian.panels
        if hybrid.lagrangian.blobs.numBlobs != 0:
            blobsFile << hybrid.lagrangian.blobs
        
        
    if timeStep % postProcStep == 0:
        
        # Export data
        eulerianFile << hybrid.multiEulerian['eulerian']
        panelFile << hybrid.lagrangian.panels
        if hybrid.lagrangian.blobs.numBlobs != 0:
            blobsFile << hybrid.lagrangian.blobs

        vxL, vyL = hybrid.lagrangian.evaluateVelocity(xL,yL)
        uLagrangian.vector()[:] = np.vstack((vxL, vyL)).T.flatten()
        
        lagrangianFieldFile << (uLagrangian,lagrangian.t)
        
        
        # Calculate parameters    
        numBlobsTime[timeStep/postProcStep] =  hybrid.lagrangian.blobs.numBlobs
        totalCirculationTime[timeStep/postProcStep] =  hybrid.lagrangian.gTotal
        totalPanelCirculationTime[timeStep/postProcStep] =  hybrid.lagrangian.panels.gTotal[0]
        totalBlobCirculationTime[timeStep/postProcStep] = hybrid.lagrangian.blobs.g.sum()
        # Calculate the body forces
        CD[timeStep/postProcStep], CL[timeStep/postProcStep] = eulerian.Forces() / (0.5*vInf[0]*vInf[0]*L)
            
        # Pressure drag
        CDpres[timeStep/postProcStep] = eulerian.PressureForces()[0] / (0.5*vInf[0]*vInf[0]*L)
        # Frictional Drag
        CDfric[timeStep/postProcStep] = eulerian.FrictionalForces()[0] / (0.5*vInf[0]*vInf[0]*L)
        
        # Store the results
        np.savez(postProcSaveDir+'postProcData.npz',
                 deltaT = hybrid.lagrangian.deltaT,
                 nTimeSteps = nTimeSteps,
                 postProcStep = postProcStep,
                 timeStep = timeStep,
                 CD = CD,
                 CL = CL,
                 CDpres = CDpres,
                 CDfric = CDfric,
                 numBlobsTime = numBlobsTime,
                 totalCirculationTime = totalCirculationTime,
                 totalBlobCirculationTime=totalBlobCirculationTime,
                 totalPanelCirculationTime=totalPanelCirculationTime,
                 endTime = time.strftime("%Y-%m-%d, %H:%M:%S", time.localtime()))        
    


        print 'Number of blobs\t\t\t: %g'  % numBlobsTime[timeStep/postProcStep]
        print 'Total circulation\t\t: %g'  % totalCirculationTime[timeStep/postProcStep]



    
    print 'Number of blobs\t\t\t: %g' % hybrid.lagrangian.blobs.numBlobs
    print 'Total circulation\t\t: %g'  % hybrid.lagrangian.gTotal  
    print '--------------------------------------------\n'

#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
# Plot results

nTimeStepsBackup = np.copy(nTimeSteps)
nTimeSteps = timeStep

# Latex Text
py.rc('text', usetex=True)
py.rc('font', family='serif')

tEval = np.arange(0,(nTimeSteps/postProcStep)+1)*hybrid.lagrangian.deltaT*postProcStep
nData = tEval.shape[0]
# Number of blobs
py.figure()
py.plot(tEval,numBlobsTime[:nData])
py.xlabel(r'$t$')
py.ylabel(r'$N_{blobs}$')
py.grid()

# Total circulation
py.figure()
py.semilogy(tEval,np.abs(totalCirculationTime[:nData]-totalCirculationTime[0]))
py.xlabel(r'$t$')
py.ylabel(r'$\Delta \Gamma$')
py.grid()

# Lift coefficient
py.figure()
py.plot(tEval,CL[:nData],'b.-',label='hybrid')
py.xlabel(r'$t$')
py.ylabel(r'$C_L$')
py.axis([0., (nTimeSteps+1)*hybrid.lagrangian.deltaT, -10,10])
py.grid()

py.plot(np.arange(0,8701)*data.f.deltaT,data.f.CL[:8701],'k',label='full Eulerian')
py.axis([0,0.87,-1,3])
py.legend()

#py.xlabel(r'$t [-]$')
#py.ylabel(r'$C_L$')
#py.axis([0., (nTimeSteps+1)*eulerian.deltaT, -5,5])
#py.grid()


# Drag ceofficeint
py.figure()
py.plot(tEval,CD[:nData],'b-',label='Total')
py.plot(tEval,CDfric[:nData],'b-.',label='Friction')
py.plot(tEval,CDpres[:nData],'b--',label='Pressure')
py.xlabel(r'$t$')
py.ylabel(r'$C_D$')
py.axis([0., (nTimeSteps+1)*hybrid.lagrangian.deltaT, 0,1])
py.legend()
py.grid()

py.figure()
py.plot(tEval,CD[:nData],'b.-',label='Hybrid')
py.plot(np.arange(0,8701)*data.f.deltaT,data.f.CD[:8701],'k',label='full Eulerian')
py.xlabel(r'$t$')
py.ylabel(r'$C_D$')
py.axis([0., (nTimeSteps+1)*hybrid.lagrangian.deltaT, 0,1])
py.legend()
py.grid()
