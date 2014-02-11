#-*- coding: utf-8 -*-
__doc__ = """

ex_05
    ****** NOTE ******* :: NO MPI POSSIBLE
    
    Convection of double Monopole
    
    
Author      :   %s
First added :   %s           
Copyright   :   %s
Licence     :   %s 
"""

__author__      = "Lento Manickathan <l.manickathan@student.tudelft.nl>"
__date__        = "2013-08-13"
__copyright__   = "Copyright (C) 2013 " + __author__
__license__     = "GNU GPL version 3 or any later version"
__doc__         %= (__author__, __date__, __copyright__ , __license__)


# Packages required               

from pHyFlow.NavierStokes import NSSolver   # main Navier-Stokes solver file

import dolfin
#import numpy as np
#import pylab as py
#import time


# Keyword for ipython greedy completion: useful for advance ipython interaction
#       %config IPCompleter.greedy = True

#===========================================================================
# Defining Parameters

# Dipole parameters
Re      = 625.0 # Reynolds number #NOTE: low-Reynolds because NSSolver does
                # not have a turbulence model implemented
                 
U       = 1.0
W       = 1.0
x1,y1   = 0.1, 1.0 # Initial core location. The core is located at (0.0,2.0), (x,y) and the wall is at x = 0.0
x2,y2   = -0.1, 1.0
r0      = 0.1   # Half-core radius

omega_e = 299.528385375226 #320.0

T       = 1.0  # Non-dimensionalized time scale

CFL     = 0.75

solverType = 'chorin'
#===========================================================================

#===========================================================================
# Calculating the parameters

tend    = T*W/U     # tend
nu      = U*W/Re    # Viscosity
Umax    = 12.0*U     # Umax

#===========================================================================

       
#===========================================================================
# Define the geometry and fluid domains

# Define the location where simulation mesh file is saved
mesh = 'geometry/dipoleConvection_exteriorBC_Re-625_vertices-41k_MPI.xml.gz'

# MeshFunction
boundaryDomains = 'geometry/dipoleConvection_exteriorBC_Re-625_vertices-41k_MPI_facet_region.xml.gz'

numVerticesK = 41
#===========================================================================

#===========================================================================
# Double-Monopole exact solution

# Define the double monopole expression
u0_doubleMonopole = dolfin.Expression(('-0.5*fabs(omega_e)*(x[1]-y1)*exp(-(pow(x[0]-x1,2)+pow(x[1]-y1,2))/pow(r0,2))\
                                        +0.5*fabs(omega_e)*(x[1]-y2)*exp(-(pow(x[0]-x2,2)+pow(x[1]-y2,2))/pow(r0,2))', # u: x-dir
                                        '0.5*fabs(omega_e)*(x[0]-x1)*exp(-(pow(x[0]-x1,2)+pow(x[1]-y1,2))/pow(r0,2))\
                                        -0.5*fabs(omega_e)*(x[0]-x2)*exp(-(pow(x[0]-x2,2)+pow(x[1]-y2,2))/pow(r0,2))'),# v: y-dir
                                        omega_e=omega_e, x1=x1, y1=y1, x2=x2, y2=y2, r0=r0) # constants
                    
#===========================================================================
                    
#===========================================================================
# Initialize the Navier-Stokes solver

# Initialize solver
nssolver = NSSolver(mesh, boundaryDomains, nu, Umax, CFL, noSlipBC=False)#, interpolationMesh_bounds=meshBounds, interpolationMesh_density=meshDensity)

# Initial Conditions
u0 = dolfin.interpolate(u0_doubleMonopole, nssolver.V) # Velocity
p0 = dolfin.Function(nssolver.Q) # Pressure field

# Assign initial conditions
nssolver.initialConditions(u0,p0)

# Standard solution
#nssolver.plot('u1','vorticity','mesh', interactive=True)
#===========================================================================


#=========================================================================== 
# Save parameters

parameterFilePath = 'results/info_dipoleConvection_Re-' + str(Re) + '_CFL-' + str(CFL) + '_vertices-' + str(numVerticesK) + 'k.txt'

saveStyle = '%15s%30s%30s\n'
with open(parameterFilePath,'wr+') as parameterFile:
    parameterFile.write(saveStyle %('Variable', 'Description', 'Data'))
    parameterFile.write(80*'-'+'\n')
    parameterFile.write(saveStyle %('Re',   'Reynolds Number', Re))
    parameterFile.write(saveStyle %('W',    'Charactersitic length', W))    
    parameterFile.write(saveStyle %('U',    'Charactersitic Velocity', U))
    parameterFile.write(saveStyle %('nu',   'Kinetmatic Viscosity', nu))
    parameterFile.write(saveStyle %('Umax', 'Maximum Velocity', Umax))
    parameterFile.write(saveStyle %('r0',   'Core radius', r0))
    parameterFile.write(saveStyle %('omega_e', 'max. vortex strength', omega_e))
    parameterFile.write(saveStyle %('x1, y1', 'Location of Core 1', (x1, y1)))
    parameterFile.write(saveStyle %('x2, y2', 'Location of Core 2', (x2, y2)))
    parameterFile.write(saveStyle %('T',    'Time Length', T))
    parameterFile.write(saveStyle %('CFL',  'Stability parameter', CFL))
    parameterFile.write(saveStyle %('solverType',  'Navier-Stokes solver', solverType))
    parameterFile.write(saveStyle %('dt',  'Time step size', nssolver.dt(0.0)))
    parameterFile.write(saveStyle %('h', 'Minimum cell size', nssolver._h))
    parameterFile.write(saveStyle %('num_vertices', 'Number of mesh vertices', nssolver.mesh.num_vertices()))
    parameterFile.write(saveStyle %('mesh', 'Mesh FilePath\t\t', mesh))
    parameterFile.write(saveStyle %('boundaryDomains', 'Mesh boundary FilePath\t\t', boundaryDomains))

parameterFile.close()
#=========================================================================== 


#===========================================================================
# Time stepping
import time

# Arbitrary number of sub-steps
nSteps = 1

# time-step size
dt = nssolver.dtMax*nSteps

# Number of iterations
iterMax = int(tend/dt)

# Boundary condition
u1_boundary = dolfin.interpolate(dolfin.Constant((0.0,0.0)), nssolver.V)

# Time stepping
for iter in range(0, iterMax):

    # Timer    
    start = time.time()
    
    # Current non-dimensional time
    T_current = (iter*dt + dt)*U/W
       
    # Solving at every step
    nssolver.step(u1_boundary, dt, nSteps)

    # Calculating simulation duration
    timeLeft = (iterMax-(iter+1))*(time.time()-start)
    timeLeftDays = timeLeft/(60.*60.*24)

    # Display the iteration parameters: current time, the nth iteration, run percentage.
    print "t = %f. \t %d of %d. \t %f%% done. Time Left: %f sec, %f days. " % (T_current, iter+1, iterMax, 100*(float(iter+1)/iterMax), timeLeft, timeLeftDays)

    #Store solution
    if iter % 10 == 0:
        nssolver.saveData('xml', (iter*nssolver.dtMax + nssolver.dtMax), 'results/dipoleConvection_Re' + str(Re) + '_CFL-' + str(CFL) + 'vertices-' + str(numVerticesK) + 'k/xmlData/')
        nssolver.saveData('pvd', (iter*nssolver.dtMax + nssolver.dtMax), 'results/dipoleConvection_Re' + str(Re) + '_CFL-' + str(CFL) + 'vertices-' + str(numVerticesK) + 'k/pvdData/') 
    
    # Plot instantanious solution
    nssolver.plot('vorticity') 
    
    
# End (and allow interaction)
print "Done."        
dolfin.interactive()
