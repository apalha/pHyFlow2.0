#!/usr/bin/env python
#-*- coding: utf-8 -*-
__doc__ = """
    ex_02
    
    DESCRIPTION:
    
        Collision of a double-Monopole (dipole), normal to the wall. The key 
        interest of this problem is dipole-Wall interaction and the generation 
        of the vorticity at the boundary layer.
        
    REFERENCES:
        Clercx, H. J. H., & Bruneau, C.-H. (2006). The normal and oblique 
        collision of a dipole with a no-slip boundary. Computers & Fluids, 
        35(3), 245–279. doi:10.1016/j.compfluid.2004.11.009
        
        
        
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
# Define the fluid and dipole parameters

# Dipole parameters
Re      = 625.0 # Reynolds number #NOTE: low-Reynolds because NSSolver does
                 # not have a turbulence model implemented
                 
U       = 1.0
W       = 1.0
x1,y1   = 0.1, 0.0  # Initial core location. Core 1 (positive)
x2,y2   = -0.1, 0.0 # Core 2 (negative)
r0      = 0.1 # Half-core radius

omega_e = 299.528385375226 #320.0

T       = 2.0  # Non-dimensionalized time scale

CFL     = 0.75

# Navier-Stokes solver algorithm
solverType = 'chorin'  # Chorin Projection scheme

#===========================================================================

#===========================================================================
# Calculating the parameters

tend    = T*W/U     # tend
nu      = U*W/Re    # Viscosity
Umax    = 12.0*U     # Umax

#===========================================================================

#===========================================================================
# Setting up the problem

# Define body geometry
mesh            = "geometry/doubleMonopole_Re-625_vertices-37k_MPI.xml.gz"               # mesh around geometry
boundaryDomains = "geometry/doubleMonopole_Re-625_vertices-37k_MPI_facet_region.xml.gz"  # mesh boundary I.Ds

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
# Initialize the problem

# Initialize the solver
nssolver =  NSSolver(mesh, boundaryDomains, nu, Umax, CFL, solverType)

# Initial conditions
u0 = dolfin.interpolate(u0_doubleMonopole, nssolver.V)
p0 = dolfin.Function(nssolver.Q)

# Initializing the initial solutions
nssolver.initialConditions(u0, p0)
nssolver.plot('u1','p1','vorticity', 'mesh', interactive=True)

#=========================================================================== 

#=========================================================================== 
# Save parameters

parameterFilePath = 'results/info_doubleMonopole_Re-' + str(Re) + '_CFL-' + str(CFL) + 'vertices-37k.txt'

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

dt = nssolver.dtMax
iterMax = int(tend/dt)

u1_boundary = dolfin.Function(nssolver.V)

for iter in range(0, iterMax):

    # Current Non-dimensional time
    T_current = (iter*nssolver.dtMax + nssolver.dtMax)*U/W    
    
    # Solving at every step
    nssolver.step(u1_boundary,dt)

    # Plot instantanious solution
    #nssolver.Plot('u1')

    # Display information    
    print "t = %f. \t %d of %d. \t %f%% done." % (T_current, iter+1, iterMax, 100*(float(iter+1)/iterMax))

    #Store solution
    nssolver.saveData('xml', (iter*nssolver.dtMax + nssolver.dtMax), 'results/ex_03/doubleMonopole_Re-' + str(Re) + '_CFL-' + str(CFL) + 'vertices-37k/xmlData/')
    
    if iter % 5 == 0:
        nssolver.saveData('pvd', (iter*nssolver.dtMax + nssolver.dtMax), 'results/ex_03/doubleMonopole_Re-' + str(Re) + '_CFL-' + str(CFL) + 'vertices-37k/') 
        
print "Done."        

# End (and allow interaction)
dolfin.interactive()

#===========================================================================
