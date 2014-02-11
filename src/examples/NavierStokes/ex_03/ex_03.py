#!/usr/bin/env python
#-*- coding: utf-8 -*-
__doc__ = """
ex_03

DESCRIPTION:
    Generate solutions of cylinder in free-stream. Simulating the flow of
    around the impulsively started cylinder subjected to uniform flow and
    artificial pertubation of rotating cylinder:
    
        between T = 3.0 and T < 3.5:
            
                u_rotation = + 0.15 * U_freestream
                
        between T = 3.5 and T < 5.0:
        
                u_rotation = -0.25 * U_freestream
                
REFERENCES:

    VALIDATION RESULTS: 
    
    ----- Vortex method ----
    Koumoutsakos, P., & Leonard, A. (1995). High-resolution simulations of the 
    flow around an impulsively started cylinder using vortex methods. Journal 
    of Fluid Mechanics, 1–38. Retrieved from 
    http://journals.cambridge.org/production/action/cjoGetFulltext?fulltextid=354808
    
    ----- Finite Different and Spectral method  ----- 
    Johnston, H., & Liu, J.-G. (2004). Accurate, stable and efficient Navier–Stokes 
    solvers based on explicit treatment of the pressure term. Journal of
    Computational Physics, 199(1), 221–259. doi:10.1016/j.jcp.2004.02.009
    
    ----- Finite Volume ------
    Chassaing, P. (1986). Numerical study and physical analysis of the pressure 
    and velocity fields in the near wake of a circular cylinder, 166.
    
    ARTIFICIAL PERTUBATION:
    
    Lecointe, Y., & Piquet, J. (1984). On the use of several compact methods 
    for the study of unsteady incompressible viscous flow round a circular 
    cylinder. Computers & fluids, 12(4), 255–280. Retrieved from 
    http://www.sciencedirect.com/science/article/pii/0045793084900094
    
Author      :   %s
First added :   %s           
Copyright   :   %s
Licence     :   %s     
"""

__author__      = "Lento Manickathan <l.manickathan@student.tudelft.nl>"
__date__        = "2013-07-22"
__copyright__   = "Copyright (C) 2013 " + __author__
__license__     = "GNU GPL version 3 or any later version"
__doc__         %= (__author__, __date__, __copyright__ , __license__)


"""
    Reviews:

"""

# Packages required to solve the navier-stokes problem
from NavierStokes import NSSolver   # main Navier-Stokes solver file
import dolfin                       # main module for FEniCS/DOLFIN packages

# Keyword for ipython greedy completion: useful for advance ipython interaction
#       %config IPCompleter.greedy = True

#===========================================================================
# Define the fluid and geometry parameters

# Define fluid parameter
Re  = 550.0    # Reynolds number #NOTE: low-Reynolds because NSSolver does
                # not have a turbulence model implemented
                
U   = 1.0       # Free-stream velocity. Uniform parallel flow in the x-direction
                # Note: y-dir is zero, so it is not defined.
               
D   = 2.0       # Characteristic Length: the diameter of the cylinder. The
                # meshing of this body should be equal
                
T   = 100.0       # Non-dimensionalized time length. This parameter determines
                # the length of the simulation.

CFL = 0.75      # Simulation CFL number. For explicit time marching scheme
                # we need CFL less than 1.0

solverType =  'chorin'  # To solve the navier-stokes problem, we would like
                        # to use the chorin's projection scheme 'chorin'. All
                        # the solver choices are found in './NavierStokes/solvers/'.
#===========================================================================


#===========================================================================
# Compute the simulation parameters

# Compute parameters
nu  = U*D/Re    # Fluid kinematic viscosity. Function of Reynolds number Re,
                # characteristic velocity U and characteristic length D
                
R   = 0.5*D     # Radius of the cylinder     
                       
Umax = 2.5*U    # Maximum fluid velocity. From potential flow around cylinder
                # in free-stream flow, the maximum fluid velocity is at the 
                # top of the cylinder and is 2 x U_freestream. The maximum
                # fluid velocity is to satisfy the CFL stability criterion and
                # to determine the time step value dt.      
                
t_final = (T*R)/U   # Simulation end time. Is a function of non-dimensionalized
                    # time length, charactersitic length D, and characterstic
                    # velocity U.               
#===========================================================================
      
       
#===========================================================================
# Define the geometry and fluid domains

# Note: The following variables carry the location of the fluid mesh and 
# boundary function files. The fluid mesh is generated using 'Gmsh' and was
# converted to 'xml.gz' using dolfin-convert. The 'mesh' is simply the mesh
# of the fluid domain around the cylinder. The 'boundaryDomains' is a mesh
# function file that carry the identification of the domain boundaries. This
# should have been pre-defined when meshing in 'Gmsh'. The boundary [2] is the
# no-slip boundary and the boundary [3] is the outer boundary of the fluid
# domain where the far-field boundary conditions are prescribed.

# Define the location where simulation mesh file is saved
mesh = "geometry/cylinder2D_Re-550_vertices-33k_MPI.xml.gz"     # mesh of the simulation fluid
                                                                 # around the geometry of interest.                             
# Define the location where domain boundary function file is saved.
boundaryDomains = "geometry/cylinder2D_Re-550_vertices-33k_MPI_facet_region.xml.gz"    # mesh boundary I.Ds

numVerticesK = 33 # thousands

#===========================================================================


#===========================================================================
# Initialize the Navier-Stokes solver

# Initialize the NSSolver. The input parameters are the mesh 'mesh', boundary mesh
# function 'boundaryDomains', viscosity 'nu', maximum velocity 'Umax' and 
# the solver scheme type 'chorin'. The NSSolver.__init__ initializes the 
# mesh functions that are needed for inputing the initial velocity, pressure
# conditions.
nssolver = NSSolver(mesh, boundaryDomains, nu, Umax, CFL, solverType)

# Inputing the initial conditions into the solver. This will assign the u0
# and the p0 inside the solver.
u0 = dolfin.Function(nssolver.V) # u,v = zero
p0 = dolfin.Function(nssolver.Q) # p = zero

nssolver.initialConditions(u0, p0)    

# Plot initial conditions
#nssolver.plot('u1', 'p1', 'vorticity', 'mesh', interactive=True)
#===========================================================================


#===========================================================================
# Save parameters to txt file

# Location of the save
parameterFilePath = 'results/ex_03/test_info_ex_03_Re-' + str(Re) + '_CFL-' + str(CFL) + '_vertices-' + str(numVerticesK) + 'k.txt'

# Save format
saveStyle = '%15s%30s%30s\n'

# Saving the data
with open(parameterFilePath,'wr+') as parameterFile:
    parameterFile.write(saveStyle %('Variable', 'Description', 'Data'))
    parameterFile.write(80*'-'+'\n')
    parameterFile.write(saveStyle %('Re',   'Reynolds Number', Re))
    parameterFile.write(saveStyle %('D',    'Charactersitic length', D))    
    parameterFile.write(saveStyle %('U',    'Charactersitic Velocity', U))
    parameterFile.write(saveStyle %('nu',   'Kinetmatic Viscosity', nu))
    parameterFile.write(saveStyle %('Umax', 'Maximum Velocity', Umax))
    parameterFile.write(saveStyle %('T',    'Time Length', T))
    parameterFile.write(saveStyle %('CFL',  'Stability parameter', CFL))
    parameterFile.write(saveStyle %('solverType',  'Navier-Stokes solver', solverType))
    parameterFile.write(saveStyle %('dt',  'Time step size', nssolver.dt(0)))
    parameterFile.write(saveStyle %('h', 'Minimum cell size', nssolver._h))
    parameterFile.write(saveStyle %('num_vertices', 'Number of mesh vertices', nssolver.mesh.num_vertices()))
    parameterFile.write(saveStyle %('mesh', 'Mesh FilePath\t\t', mesh))
    parameterFile.write(saveStyle %('boundaryDomains', 'Mesh boundary FilePath\t\t', boundaryDomains))

parameterFile.close() # Close the file
#===========================================================================


#===========================================================================
# Iterative time stepping with NSSolver

import time

# Determine the time step size. 
dt = nssolver.dtMax # For now, we wish to used the nssolver's maximum
                    # time step size 'dtMax'

# Determine the maximum number of time steps. 
iterMax = int(t_final/dt) # Number of iteration is determined by the number of files

# Iterating
for iter in range(0, iterMax):
        
    # Timer        
    start = time.time()
    
    # Current time
    T_current = (iter*dt + dt)*U/(D*0.5)
    
    # Stepping
    nssolver._step_standard(U, T=T_current)


    # Timer - calculating duration of simulation
    timeLeft = (iterMax-(iter+1))*(time.time()-start)
    timeLeftDays = timeLeft/(60.*60.*24)

    # Display the iteration parameters: current time, the nth iteration, run percentage.
    print "t = %f. \t %d of %d. \t %f%% done. Time Left: %f sec, %f days. " % (T_current, iter+1, iterMax, 100*(float(iter+1)/iterMax), timeLeft, timeLeftDays)

    # Store pressure and velocity data
    #nssolver.saveData('xml', T_current, 'results/ex_03/ex_03_Re-' + str(Re) + '_CFL-' + str(CFL) + 'vertices-' + str(numVerticesK) + 'k/xmlData/')

    # Save the Plot data for further plotting in advanced plotting application
    # such as 'PARAVIEW'. For such evaluations only plot of every 50-100
    # iterations is needed.
    if iter % 10 == 0:
        # Store solution at every time-step.
        nssolver.saveData('pvd', T_current, 'results/ex_03/ex_03_Re-' + str(Re) + '_CFL-' + str(CFL) + 'vertices-' + str(numVerticesK) + 'k/pvdData/')

    # Plot instantanious solution
    nssolver.plot('vorticity', 'p1')
    
# END
dolfin.interactive()
#===========================================================================

