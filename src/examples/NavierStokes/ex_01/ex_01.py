#!/usr/bin/env python
#-*- coding: utf-8 -*-
__doc__ = """
    ex_01
    
    DESCRIPTION:
        Example 01 showing how the Navier-Stokes solver works. The Navier-Stokes
        solver for time-stepping with just the Dirichlet boundary conditions
        is named "NSSolver.py" and is located in the folder 'NavierStokes'.
        
        Example 01 involves simulating a cylinder submerged in a laminar 
        freestream flow. To simulate such conditions, one need to prescribe 
        a no-slip Dirichlet B.C at the cylinder wall [boundary Domain (2)]
        and a far-field free-steam flow condition at the edge of the 
        fluid domain [boundary Domain [3]]. The dirichlet B.C is obtained
        from the exact solutions generated before and stored in './results/'.
        
        In order to solver the incompressible Navier-Stokes equation,
        the Finite-Element method uses Chorin's projection method to find the
        solution of the next time step.
        

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
import glob                         # filename globbing utility

# Keyword for ipython greedy completion: useful for advance ipython interaction
#       %config IPCompleter.greedy = True

#===========================================================================
# Define the fluid and geometry parameters

# Define fluid parameter
Re  = 100.0     # Reynolds number #NOTE: low-Reynolds because NSSolver does
                # not have a turbulence model implemented
                
U   = 1.0       # Free-stream velocity. Uniform parallel flow in the x-direction
                # Note: y-dir is zero, so it is not defined.
               
D   = 0.1       # Characteristic Length: the diameter of the cylinder. The
                # meshing of this body should be equal
                
T   = 100.0     # Non-dimensionalized time length. This parameter determines
                # the length of the simulation. For a cylinder in uniform
                # parallel flow, Karmen-Vortex street starts appearing around
                # T = 80

CFL = 0.25      # Simulation CFL number. For explicit time marching scheme
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

tend    = T*(D/2)/U     # Simulation end time. Is a function of non-dimensionalized
                        # time length, charactersitic length D, and characterstic
                        # velocity U.
                        
Umax = 2.5*U    # Maximum fluid velocity. From potential flow around cylinder
                # in free-stream flow, the maximum fluid velocity is at the 
                # top of the cylinder and is 2 x U_freestream. The maximum
                # fluid velocity is to satisfy the CFL stability criterion and
                # to determine the time step value dt.      
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
mesh = "geometry/cylinder2D_fine.xml.gz"     # mesh of the simulation fluid
                                             # around the geometry of interest.
                                                            
# Define the location where domain boundary function file is saved.
boundaryDomains = "geometry/cylinder2D_fine_facet_region.xml.gz"    # mesh boundary I.Ds
#===========================================================================

#===========================================================================
# Load the saved exact solutions

# Location of the velocity and pressure data .xml.gz files
dataPath_u  = "results/ex_01/data_u_t*.xml.gz" # velocity
dataPath_p  = "results/ex_01/data_p_t*.xml.gz" # pressure

# Load the data files into a list and sort it alphanumerically
# Note:  'glob.glob' return list of paths matching the 'dataPath' pattern
dataFileList_u  = sorted(glob.glob(dataPath_u)) # velocity
dataFileList_p  = sorted(glob.glob(dataPath_p)) # pressure

# Load the mesh used to solve the exact solution
# Note: this is needed to load the solution to 'DOLFIN' function format
mesh_exact = dolfin.Mesh('geometry/cylinder2D_fine.xml.gz')

# Plot the mesh of the exact solution
# Note: key is the plot reference I.D. In order to replot in the same plot
# same key must be used. Also, in order to interact with the plot, interactive
# has to be set to True.
dolfin.plot(mesh_exact, title='exact_mesh', key='exact_mesh', interactive=True)


# Define the function spaces of the exact solution
V_exact     = dolfin.VectorFunctionSpace(mesh_exact, "CG", 2) # velocity function space
Q_exact     = dolfin.FunctionSpace(mesh_exact, "CG", 1)       # pressure function space  
#===========================================================================

#===========================================================================
# Initialize the Navier-Stokes solver

# Initialize the NSSolver. The input parameters are the mesh 'mesh', boundary mesh
# function 'boundaryDomains', viscosity 'nu', maximum velocity 'Umax' and 
# the solver scheme type 'chorin'. The NSSolver.__init__ initializes the 
# mesh functions that are needed for inputing the initial velocity, pressure
# conditions.
nssolver = NSSolver(mesh, boundaryDomains, nu, Umax, CFL, solverType)

# Initial conditions: intial velocity u0 and initial pressure p0
# Load the initial solution using dataFileList and function space V_exact, Q_exact
u0_exact = dolfin.Function(V_exact, dataFileList_u[0]) # ux, uy
p0_exact = dolfin.Function(Q_exact, dataFileList_p[0]) # p

# Intepolate the exact solution to the function spaces V (vector) and Q (scalar)
# of the solver mesh.
u0_interpolated = dolfin.interpolate(u0_exact, nssolver.V)
p0_interpolated = dolfin.interpolate(p0_exact, nssolver.Q)

# Inputing the initial conditions into the solver. This will assign the u0
# and the p0 inside the solver.
nssolver.initialConditions(u0_interpolated, p0_interpolated)
#===========================================================================

#===========================================================================
# Plotting the variables

# Plot the initial conditions. NSSolver has a 'Plot' function where it is able
# to plot any plottable parameters such as u0, u1, p0, p1 and the mesh. Plots
# the variables stored in nssolver with the proper argument name.
nssolver.plot('mesh','u0','p0', interactive=True) # input any number of arguments with the proper variables
#===========================================================================

#===========================================================================
# Single time stepping

# Determine the time step size. For now, we wish to used the nssolver's maximum
# time step size 'dtMax'
dt = nssolver.dtMax

# Load the boundary conditions from the exact solution of the next step
# The exact solutions were generated before hand, and only the b.c
# at the outer domain is need to calculated the fluid parameters in the
# fluid domain
u1_exact = dolfin.Function(V_exact, dataFileList_u[1])
p1_exact = dolfin.Function(Q_exact, dataFileList_p[1])

# Interpolating the solution to the solver function space.
u1_interpolated = dolfin.interpolate(u1_exact, nssolver.V)
p1_interpolated = dolfin.interpolate(p1_exact, nssolver.Q)

# Feeding the velocity solution into the solver. The solver will extract
# the dirichlet velocity b.c where the solver's mesh boundary is.    
nssolver.step(u1_interpolated, dt)

# Plot the solutions of the new time step, velocity, pressure and vorticity
nssolver.plot('u1','p1','vorticity', interactive=True)
#===========================================================================


#===========================================================================
# Storing the solutions

# .xml files for post-processing and post-calculations
# Store solution at every time-step. This save solution can be later used
# for further calculations. The 'Save_XMLDATA' saves the velocity and 
# pressure solution of the new time step u1, p1 in '.xml.gz' compressed
# format to the any location.
#       nssolver.Save_XMLData(iter*nssolver.dtMax,'./results/xmlData/')

# .vtu/.pvd files for plotting
# Save the Plot data for further plotting in advanced plotting application
# such as 'PARAVIEW'. For such evaluations only plot of every 50-100
# iterations is needed.
#           nssolver.Save_Plot('./results/')
#===========================================================================

#===========================================================================
# Iterative time stepping with NSSolver

# Reset the initial conditions
# Initial conditions: intial velocity u0 and initial pressure p0
# Load the initial solution using dataFileList and function space V_exact, Q_exact
u0_exact = dolfin.Function(V_exact, dataFileList_u[0]) # ux, uy
p0_exact = dolfin.Function(Q_exact, dataFileList_p[0]) # p
# Intepolate the exact solution to the function spaces V (vector) and Q (scalar)
# of the solver mesh.
u0_interpolated = dolfin.interpolate(u0_exact, nssolver.V)
p0_interpolated = dolfin.interpolate(p0_exact, nssolver.Q)
# Inputing the initial conditions into the solver. This will assign the u0
# and the p0 inside the solver.
nssolver.initialConditions(u0_interpolated, p0_interpolated)

# Determine the time step size. For now, we wish to used the nssolver's maximum
# time step size 'dtMax'
dt = nssolver.dtMax

# Determine the maximum number of time steps. Number of iteration is
# determined by the number of files
iterMax = len(dataFileList_u) - 1

# Iterating through all the time steps. From 1 to maximum iterations.
for iter in range(1, iterMax):
    
    # Load the boundary conditions from the exact solution of the next step
    # The exact solutions were generated before hand, and only the b.c
    # at the outer domain is need to calculated the fluid parameters in the
    # fluid domain
    u1_exact = dolfin.Function(V_exact, dataFileList_u[iter])
    p1_exact = dolfin.Function(Q_exact, dataFileList_p[iter])
    
    # Interpolating the solution to the solver function space.
    u1_interpolated = dolfin.interpolate(u1_exact, nssolver.V)
    p1_interpolated = dolfin.interpolate(p1_exact, nssolver.Q)
    
    # Feeding the velocity solution into the solver. The solver will extract
    # the dirichlet velocity b.c where the solver's mesh boundary is.    
    nssolver.step(u1_interpolated, dt)

    # Plot the solutions of the new time step, velocity, pressure and vorticity
    nssolver.plot('u1','p1','vorticity') # Note: interactive=True only after the simulation run
        
    # Display the iteratio parameters: current time, the nth iteration, run percentage.
    print "t = %f. \t %d of %d. \t %f%% done." % (iter*nssolver.dtMax, iter, iterMax-1, 100*(float(iter)/(iterMax-1)))

    # Store solution at every time-step. This save solution can be later used
    # for further calculations. The 'Save_XMLDATA' saves the velocity and 
    # pressure solution of the new time step u1, p1 in '.xml.gz' compressed
    # format
    #       nssolver.Save_XMLData(iter*nssolver.dtMax)
    
    # Save the Plot data for further plotting in advanced plotting application
    # such as 'PARAVIEW'. For such evaluations only plot of every 50-100
    # iterations is needed.
    #       if iter % 50 == 0:
    #           nssolver.save_solution()

# Allow Plot interaction after the end of simuluation.
dolfin.interactive()
#===========================================================================

#=========================================================================== 
# Calculate the error between exact solution

# Caculated the error
uError = u1_exact -  nssolver.u1    # velocity
pError = p1_exact - nssolver.p1     # pressure

# Plot the error
dolfin.plot(uError, title='error in velocity') # Error in velocity
dolfin.plot(pError, title='error in pressure') # Error in pressure

# Calculate the relative Error
uRelError = abs(u1_exact - nssolver.u1)/abs(u1_exact.vector().max())
pRelError = abs(p1_exact - nssolver.p1)/abs(p1_exact.vector().max())

# Plot the error
dolfin.plot(uRelError, title='Relative error in velocity') # Error in velocity
dolfin.plot(pRelError, title='Relative error in pressure') # Error in pressure

# Interactive on
dolfin.interactive()
#=========================================================================== 
