__doc__ = """

ex_05
    ****** NOTE ******* :: NO MPI POSSIBLE
    
    Coupling Navier-Stokes to Vortex method (Just one-way). Time-stepping of 
    the navier-stokes with the Dirichlet B.C of the vortex method.
    
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
from NavierStokes import NSSolver   # main Navier-Stokes solver file
import dolfin
import numpy as np

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

CFL     = 0.95

solverType = 'chorin'
#===========================================================================

#===========================================================================
# Calculating the parameters

tend    = T*W/U     # tend
nu      = U*W/Re    # Viscosity
Umax    = 11.0*U     # Umax

#===========================================================================

#===========================================================================
# Define the geometry and fluid domains

# Mesh
mesh = 'geometry/dipoleConvection_coupled_Re-625_vertices-6k_MPI.xml.gz'
# MeshFunction
boundaryDomains = 'geometry/dipoleConvection_coupled_Re-625_vertices-6k_MPI_facet_region.xml.gz'

numVerticesK = 6 # number of vertices
#===========================================================================

#===========================================================================
# Initialize the Navier-Stokes solver

# Initialize solver
nssolver = NSSolver(mesh, boundaryDomains, nu, Umax, CFL, noSlipBC=False)#, interpolationMesh_bounds=meshBounds, interpolationMesh_density=meshDensity)

# Initial Conditions
u0 = dolfin.Function(nssolver.V) # Velocity
p0 = dolfin.Function(nssolver.Q) # Pressure field

# Assign initial conditions
nssolver.initialConditions(u0,p0)

# Standard solution
#nssolver.plot('u1','vorticity','mesh', interactive=True)
#===========================================================================

#===========================================================================
# Export domain boundary coordinates - [for Vortex Method]

xboundary,yboundary = nssolver.boundary_DOFCoordinates
# boundary coordinates to txt file
#np.savetxt('results/ex_06/NS_boundaryCoordinates_Re-' + str(Re) + '_vertices' + str(numVerticesK) + 'k.txt', nssolver.boundary_DOFCoordinates)

#===========================================================================


#===========================================================================
# Load dirichlet b.c data

import scipy.io as sio

#import glob

# Nblobs = 400 x 3*400
# File location
#dirichletBC_nBlobs400_dataLocation = '/home/lento/Documents/programs/vortexMethods/' +\
#                                     'vortexlib2d_matlab/testCases/doubleMonopole/'  +\
#                                     'generate_NS_DirichletBC_results/old/data_dirichletBC.mat'

# Loading the .mat file
#dirichletBC_nBlobs400_data = sio.loadmat(dirichletBC_nBlobs400_dataLocation)

# Extracting the data and storing to variables
#dirichletBC_nBlobs400_T  = dirichletBC_nBlobs400_data[dirichletBC_nBlobs400_data.keys()[0]]['T'][0][0][0]
#dirichletBC_nBlobs400_vx = dirichletBC_nBlobs400_data[dirichletBC_nBlobs400_data.keys()[0]]['vx'][0][0]
#dirichletBC_nBlobs400_vy = dirichletBC_nBlobs400_data[dirichletBC_nBlobs400_data.keys()[0]]['vy'][0][0]

# Nblobs = 1000
nBlobs = 1000
# Discrete vorticity field: Directly from blob circulation
wDiscrete_nBlobs1000_dataLocation = '/home/lento/Documents/programs/vortexMethods/' +\
                                    'vortexlib2d_matlab/testCases/doubleMonopole' +\
                                    '/generate_NS_DirichletBC_results/data_nBlobs-1000_wDiscrete_T-*.mat'

# vx induction on the domain boundary
vx_nBlobs1000_dataLocation = '/home/lento/Documents/programs/vortexMethods/' +\
                             'vortexlib2d_matlab/testCases/doubleMonopole/' +\
                             'generate_NS_DirichletBC_results/data_nBlobs-1000_vx.txt'

# vy induction on the domain boundary
vy_nBlobs1000_dataLocation = '/home/lento/Documents/programs/vortexMethods/' +\
                             'vortexlib2d_matlab/testCases/doubleMonopole/' +\
                             'generate_NS_DirichletBC_results/data_nBlobs-1000_vy.txt'

                             
T_nBlobs1000_dataLocation = '/home/lento/Documents/programs/vortexMethods/' +\
                             'vortexlib2d_matlab/testCases/doubleMonopole/' +\
                             'generate_NS_DirichletBC_results/data_nBlobs-1000_T.txt'

# Load the data                             
vx_nBlobs1000 = np.loadtxt(vx_nBlobs1000_dataLocation) # vx
vy_nBlobs1000 = np.loadtxt(vy_nBlobs1000_dataLocation) # vx
T_nBlobs1000  = np.loadtxt(T_nBlobs1000_dataLocation) # T - non-dimensional

#===========================================================================


#===========================================================================
# Time - Stepping

import time

# Dt of vortex method - non-dimensional to dimensional
#dt = (dirichletBC_T[1] - dirichletBC_T[0])*W/U  # nBlobs = 400
dt = (T_nBlobs1000[1] - T_nBlobs1000[0])*W/U # nBlobs = 1000

# Number of sub-steps
nSteps = int(np.ceil(dt/nssolver.dtMax))

# Number of total iterations
iterMax = int(tend/dt)

for iter in range(0,iterMax-1):

    # Timer    
    start = time.time()
    
    # Current non-dimensional time
    T_current = (iter*dt + dt)*U/W
    
    # Package boundary condition
    #u1_boundary = np.array((dirichletBC_vx[iter,:],dirichletBC_vy[iter,:])) # nBlobs = 400
    u1_boundary = np.array((vx_nBlobs1000[iter,:], vy_nBlobs1000[iter,:])) # nBlobs = 1000

    # Solving at every step
    #nssolver.step(u1_boundary, (dirichletBC_T[iter+1] - dirichletBC_T[iter]), nSteps)  # nBlobs = 400
    nssolver.step(u1_boundary, (T_nBlobs1000[iter+1] - T_nBlobs1000[iter]), nSteps)  # nBlobs = 1000
    
    # Calculating simulation duration
    timeLeft = (iterMax-(iter+1))*(time.time()-start)
    timeLeftDays = timeLeft/(60.*60.*24)

    # Display the iteration parameters: current time, the nth iteration, run percentage.
    print "t = %f. \t %d of %d. \t %f%% done. Time Left: %f mins, %f hrs, %f days. " % (T_current, iter+1, iterMax, 100*(float(iter+1)/iterMax), timeLeft/(60.), timeLeft/(60.*60.), timeLeft/(60.*60.*24))

    # Plotting
    nssolver.plot('u1','p1','vorticity')
    
    #Store solution
    #nssolver.saveData('xml', dirichletBC_T[iter+1], 'results/ex_06/dipoleConvection_Re' + str(Re) + '_CFL-' + str(CFL) + 'vertices-' + str(numVerticesK) + 'k/xmlData/')
    #nssolver.saveData('pvd', dirichletBC_T[iter+1], 'results/ex_06/dipoleConvection_Re' + str(Re) + '_CFL-' + str(CFL) + 'vertices-' + str(numVerticesK) + 'k/pvdData/') 
   
    nssolver.saveData('xml', T_nBlobs1000[iter+1], 'results/ex_06/dipoleConvection_nBlobs-' + str(nBlobs) + '_Re-' + str(Re) + '_CFL-' + str(CFL) + 'vertices-' + str(numVerticesK) + 'k/xmlData/')
    nssolver.saveData('pvd', T_nBlobs1000[iter+1], 'results/ex_06/dipoleConvection_nBlobs-' + str(nBlobs) + '_Re-' + str(Re) + '_CFL-' + str(CFL) + 'vertices-' + str(numVerticesK) + 'k/pvdData/') 
      
   
    # Modifiy view angle    
    if iter == 0:
        dolfin.interactive()
    
# DONE.
print "Done."
dolfin.interactive()
#===========================================================================


#=========================================================================== 
# Save parameters

parameterFilePath = 'results/ex_05/info_dipoleConvection_Re-' + str(Re) + '_CFL-' + str(CFL) + '_vertices-' + str(numVerticesK) + 'k.txt'

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