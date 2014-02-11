#!/usr/bin/env python
#-*- coding: utf-8 -*-
__doc__ = """
ex_04

DESCRIPTION:
    Post-Processing ex_03. Calculating the lift and the drag of the cylinder in
    free-stream.

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
#from NavierStokes import NSSolver   # main Navier-Stokes solver file
#import dolfin                   
import numpy as np
import os


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
                
#T   = 1.0#np.array([5.0,6.0]) # Non-dimensionalized time length. Coarse and Fine simulation parameters

#===========================================================================


#===========================================================================
# Compute the simulation parameters

# Compute parameters
nu  = U*D/Re    # Fluid kinematic viscosity. Function of Reynolds number Re,
                # characteristic velocity U and characteristic length D
                
R   = 0.5*D                

#t_final = T*(D/2)/U     # Simulation end time. Is a function of non-dimensionalized
                        # time length, charactersitic length D, and characterstic
                        # velocity U.
                        
Umax = 2.5*U    # Maximum fluid velocity. From potential flow around cylinder
                # in free-stream flow, the maximum fluid velocity is at the 
                # top of the cylinder and is 2 x U_freestream. The maximum
                # fluid velocity is to satisfy the CFL stability criterion and
                # to determine the time step value dt.    
                
                

#===========================================================================
      

#===========================================================================
# Simulation Run parameters

# Time step - dimensional
#dt = [0.000652520809173, 0.000501742083944] 
#
## Mesh file
#meshFile = ['geometry/cylinder2D_Re-1000_vertices-26k_MPI.xml.gz',
#            'geometry/cylinder2D_smallFineDomain_Re-1000_vertices-195k_MPI.xml.gz']
#           
## Boundary mesh file           
#boundaryDomainsFile = ['geometry/cylinder2D_Re-1000_vertices-26k_MPI_facet_region.xml.gz',
#                       'geometry/cylinder2D_smallFineDomain_Re-1000_vertices-195k_MPI_facet_region.xml.gz']  
#
## Velocity XML data path
#dataPath_u = ['results/ex_03/ex_03_Re-1000.0_CFL-0.1vertices-26k/xmlData/data_u_t*.xml.gz',
#              'results/ex_03/ex_03_Re-1000.0_CFL-0.1vertices-195k/xmlData/data_u_t*.xml.gz']
#
## Pressure XML data path              
#dataPath_p = ['results/ex_03/ex_03_Re-1000.0_CFL-0.1vertices-26k/xmlData/data_p_t*.xml.gz',
#              'results/ex_03/ex_03_Re-1000.0_CFL-0.1vertices-195k/xmlData/data_p_t*.xml.gz']

dt = np.array([0.00185345000003]) #np.array([0.00180774474715, 0.000348560617885])

saveFrequency = 1 #[10,1]

readFrequency = 100

# Mesh file
meshFile = ["geometry/cylinder2D_Re-550_vertices-33k_MPI.xml.gz"]
           
# Boundary mesh file           
boundaryDomainsFile = ["geometry/cylinder2D_Re-550_vertices-33k_MPI_facet_region.xml.gz"]

# Velocity XML data path
dataPath_u = ['results/ex_03/ex_03_Re-550.0_CFL-0.75vertices-33k/xmlData/data_u_t*.xml.gz']
              

# Pressure XML data path              
dataPath_p = ['results/ex_03/ex_03_Re-550.0_CFL-0.75vertices-33k/xmlData/data_p_t*.xml.gz']

## Mesh file
#meshFile = ["geometry/cylinder2D_Re-550_vertices-36k_MPI.xml.gz",
#            "geometry/cylinder2D_Re-550_vertices-154k_MPI.xml.gz"]
#           
## Boundary mesh file           
#boundaryDomainsFile = ["geometry/cylinder2D_Re-550_vertices-36k_MPI_facet_region.xml.gz",
#                       "geometry/cylinder2D_Re-550_vertices-154k_MPI_facet_region.xml.gz"]
#
## Velocity XML data path
#dataPath_u = ['results/ex_03/ex_03_Re-550.0_CFL-0.75vertices-36k/xmlData/data_u_t*.xml.gz',
#              'results/ex_03/ex_03_Re-550.0_CFL-0.75vertices-154k/xmlData/data_u_t*.xml.gz']
#              
#
## Pressure XML data path              
#dataPath_p = ['results/ex_03/ex_03_Re-550.0_CFL-0.75vertices-36k/xmlData/data_p_t*.xml.gz',
#              'results/ex_03/ex_03_Re-550.0_CFL-0.75vertices-154k/xmlData/data_p_t*.xml.gz']
                       
#===========================================================================                       
                       
                       
#===========================================================================
# Drag Calculation Function

import glob     # filename globbing utility
import dolfin   # main module for FEniCS/DOLFIN packages

# Kinetic energy of the pressure term
def epsilon(u):
    "Friction term"
    return 0.5*(dolfin.grad(u) + (dolfin.grad(u)).T) 

def calc_Forces(runNumber):
    """ 
        Calculate the drag
        
        Input:
            runNumber : [0 or 1] ::     0 = coarse run
                                        1 = fine run    
    """
    
    # Display info
    print "\nRun number: %d\n" % runNumber
    
    ### Pre-loading
    
    # Load mesh, boundary function
    #mesh = dolfin.Mesh(meshFile[runNumber]) # mesh
    #boundaryDomains = dolfin.MeshFunction('size_t', mesh,boundaryDomainsFile[runNumber]) # boundaries
    mesh = dolfin.Mesh(meshFile[runNumber]) # mesh
    boundaryDomains = dolfin.MeshFunction('size_t', mesh,boundaryDomainsFile[runNumber]) # boundaries
    
    # Calculate domain parameters
    n   = dolfin.FacetNormal(mesh) # normal vector
    ds  = dolfin.Measure("ds")[boundaryDomains] # line integrator
    
    eX = dolfin.Constant((1.0,0.0)) # x-dir unit vector     
    eY = dolfin.Constant((0.0,1.0)) # y-dir unit vector
    
    # Function Spaces
    Q = dolfin.FunctionSpace(mesh, 'CG', 1) # scalar function space: pressure, vorticity
    V = dolfin.VectorFunctionSpace(mesh, 'CG', 2) # vector function space: velocity
    
    # Load the pressure,velocity data path lists
    #dataFileList_u = sorted(glob.glob(dataPath_u[runNumber])) # velocity
    #dataFileList_p = sorted(glob.glob(dataPath_p[runNumber])) # Pressure
    dataFileList_u = sorted(glob.glob(dataPath_u[runNumber]), key=os.path.getmtime) # velocity
    dataFileList_p = sorted(glob.glob(dataPath_p[runNumber]), key=os.path.getmtime) # Pressure
    
    ###  Calculate the drag
    
    # Pre-Allocate memory
    #Cd_pressure = np.zeros(len(dataFileList_u)) # simple pressure equation
    
    # Navier-Stokes stress tensor term
    Cd_stressTensor_friction = np.zeros(len(dataFileList_u)) # friction
    Cd_stressTensor_pressure = np.zeros(len(dataFileList_u)) # pressure
    Cd_stressTensor_total    = np.zeros(len(dataFileList_u)) # total = friction + pressure
    
    Cl = np.zeros(len(dataFileList_u))
    
    # Iterating through all the data and calculating drag
    for iter in range(0,len(dataFileList_u),readFrequency):
        
        # Display info
        print "%d of %d." % (iter+1, len(dataFileList_u))
        
        # Load Data
        p = dolfin.Function(Q, dataFileList_p[iter])
        u = dolfin.Function(V, dataFileList_u[iter])
        # Calculate drag coefficient
        
        # Simple pressure equation
        #Cd_pressure[iter] = ((dolfin.assemble(dolfin.dot(p,n[0])*ds(2)))*2.0)/((U**2)*D)
        
        # Navier-Stokes stress tensor term
        Cd_stressTensor_friction[iter] = (-(dolfin.assemble(dolfin.dot(dolfin.dot( (2.0*nu*epsilon(u)) ,n),eX)*ds(2)))*2.0)/((U**2)*D)
        Cd_stressTensor_pressure[iter] = (-(dolfin.assemble(dolfin.dot(dolfin.dot( (-p*dolfin.Identity(u.cell().d)) ,n),eX)*ds(2)))*2.0)/((U**2)*D)
        Cd_stressTensor_total[iter]    = (-(dolfin.assemble(dolfin.dot(dolfin.dot( (2.0*nu*epsilon(u) - p*dolfin.Identity(u.cell().d)) ,n),eX)*ds(2)))*2.0)/((U**2)*D)

        # Lift
        Cl[iter] = (-(dolfin.assemble(dolfin.dot(dolfin.dot( (2.0*nu*epsilon(u) - p*dolfin.Identity(u.cell().d)) ,n),eY)*ds(2)))*2.0)/((U**2)*D)
        
    print 'Iteration Done.'
       # filename globbing utility
       # Cd_pressure, 
    return Cd_stressTensor_friction, Cd_stressTensor_pressure, Cd_stressTensor_total, Cl


#===========================================================================




#===========================================================================
# Main calling section

import pylab as py                         

# Calculate the drag coefficients
#coarse_Cd_pressure, coarse_Cd_stressTensor_friction, coarse_Cd_stressTensor_pressure, coarse_Cd_stressTensor_total = calc_Drag(0)
#fine_Cd_pressure, fine_Cd_stressTensor_friction, fine_Cd_stressTensor_pressure, fine_Cd_stressTensor_total = calc_Drag(1)
#
## Calculate the time length - Non-dimensionalized
#coarse_t    = (np.arange(0,len(coarse_Cd_pressure)*5)*((dt[0]*U)/R))[::5]
#fine_t      = np.arange(0,len(fine_Cd_pressure))*((dt[1]*U)/R)
## Plot the data
##t = np.arange(0,len(dataFileList_u)*5)*((dt*U)/R)
#
#py.ion()
#
## Figure 1: Simple pressure equation
#py.figure(1)
#py.plot(coarse_t, coarse_Cd_pressure,'b-', label='coarse')
#py.plot(fine_t, fine_Cd_pressure,'g-', label='fine')
#py.axis([0,coarse_t[-1],0,2])
#py.xlabel('T [-]')
#py.ylabel('Drag Coefficient $C_d$ [-]')
#py.title('From simple pressure expression')
#py.grid()
#py.legend()
#
## Figure 1: Simple pressure equation
#py.figure(2)
#py.plot(coarse_t, coarse_Cd_stressTensor_friction,'b--', label='coarse [friction]')
#py.plot(coarse_t, coarse_Cd_stressTensor_pressure,'b-.', label='coarse [pressure]')
#py.plot(coarse_t, coarse_Cd_stressTensor_total,'b-', label='coarse [total]')
#
#py.plot(fine_t, fine_Cd_stressTensor_friction,'g--', label='fine [friction]')
#py.plot(fine_t, fine_Cd_stressTensor_pressure,'g-.', label='fine [pressure]')
#py.plot(fine_t, fine_Cd_stressTensor_total,'g-', label='fine [total]')
#py.axis([0,coarse_t[-1],0,2])
#py.xlabel('T [-]')
#py.ylabel('Drag Coefficient $C_d$ [-]')
#py.title('From Navier-Stokes stress tensor term $\sigma$')
#py.grid()
#py.legend()


 
 
## Coarse
#Cd_stressTensor_friction0, Cd_stressTensor_pressure0, Cd_stressTensor_total0 = calc_Drag(0)
#t0 = (np.arange(0,len(Cd_stressTensor_friction0))*saveFrequency[0]*((dt[0]*U)/R))
## Fine
#
#Cd_stressTensor_friction1, Cd_stressTensor_pressure1, Cd_stressTensor_total1 = calc_Drag(1)
#t1 = (np.arange(0,len(Cd_stressTensor_friction1))*saveFrequency[1]*((dt[1]*U)/R))
#
#py.ion()
#py.figure(1)
#
#py.plot(t0, Cd_stressTensor_total0,'k-', label='Total - [34k vertices]')
#py.plot(t0, Cd_stressTensor_friction0,'k--', label='Friction - [34k vertices]')
#py.plot(t0, Cd_stressTensor_pressure0,'k-.', label='Pressure - [34k vertices]')
#
#py.plot(t1[::100], Cd_stressTensor_total1[::100],'b-', label='Total - [154k vertices]')
#py.plot(t1[::100], Cd_stressTensor_friction1[::100],'b--', label='Friction - [154k vertices]')
#py.plot(t1[::100], Cd_stressTensor_pressure1[::100],'b-.', label='Pressure - [154k vertices]')
#
#py.axis([0,t0[-1],0,2])
#py.xlabel('T [-]')
#py.ylabel('Drag Coefficient $C_d$ [-]')
#py.title('From Navier-Stokes stress tensor term $\sigma$')
#py.yticks(np.arange(0,2.25,0.25))
#py.grid()
#py.legend(loc=0)


Cd_stressTensor_friction, Cd_stressTensor_pressure, Cd_stressTensor_total, Cl = calc_Forces(0)
t = (np.arange(0,len(Cd_stressTensor_friction))*saveFrequency*((dt[0]*U)/R))[::readFrequency]

py.ion()
py.figure(1)

py.plot(t, Cd_stressTensor_total[::readFrequency],'k-', label='Total')
py.plot(t, Cd_stressTensor_friction[::readFrequency],'k--', label='Friction')
py.plot(t, Cd_stressTensor_pressure[::readFrequency],'k-.', label='Pressure')

py.axis([0,t[-1],0,2])
py.xlabel('T [-]')
py.ylabel('Drag Coefficient $C_d$ [-]')
py.title('Drag Coefficient')
py.xticks(np.arange(0,110,10))
py.yticks(np.arange(0,2.25,0.25))
py.grid()
py.legend(loc=0)

py.figure(2)
py.plot(t, Cl[::readFrequency],'k-')

#py.axis([0,t[-1],0,2])
py.xlabel('T [-]')
py.ylabel('Lift Coefficient $C_l$ [-]')
py.title('Lift Coefficient')
py.xticks(np.arange(0,110,10))
py.yticks(np.arange(-2.0,2.0,0.25))
py.grid()
#py.legend(loc=0)


#===========================================================================





                       
#===========================================================================
# Load data

## Load mesh
#mesh = dolfin.Mesh("geometry/cylinder2D_Re-1000_vertices-26k_MPI.xml.gz")     # mesh of the simulation fluid                                                           
#boundaryDomains = dolfin.MeshFunction('size_t', mesh, "geometry/cylinder2D_Re-1000_vertices-26k_MPI_facet_region.xml.gz") # mesh boundary I.Ds
#
#n = dolfin.FacetNormal(mesh) # normal vector
#
#boundaries = dolfin.FacetFunction("size_t", mesh)
#boundaries.set_all(0)
#boundaries.set_values(boundaryDomains.array())
#ds = dolfin.Measure("ds")[boundaries]
#
## Function Space
#Q = dolfin.FunctionSpace(mesh, 'CG', 1)
#V = dolfin.VectorFunctionSpace(mesh, 'CG', 2)
#
##X = dolfin.VectorFunctionSpace(mesh, 'CG', 1)
#
## Location of the body
##class BodyRegion(dolfin.SubDomain):
##    def inside(self, x, on_boundary):
##        return False
###        return x[0] > -0.2 and x[0] < 0.2 and \
###               x[1] > 0.0 and x[1] < 1.3
##               #x[1] > -0.2 and x[1] < 0.2
##        #return x[0] > -1.0+dolfin.DOLFIN_EPS and x[0] < 3.0-dolfin.DOLFIN_EPS and \
##        #       x[1] > -2.0+dolfin.DOLFIN_EPS and x[1] < 2.0-dolfin.DOLFIN_EPS
##                       
#               
#def epsilon(u):
#    "Friction term"
#    return 0.5*(dolfin.grad(u) + (dolfin.grad(u)).T)   
#
#eX = dolfin.Constant((1.0,0.0))             
#eY = dolfin.Constant((0.0,1.0))
#
## Set up the functions
##bodyRegion = BodyRegion()
#
## Load pressure data
#dataPath_u = "results/ex_03/ex_03_Re-1000.0_CFL-0.1vertices-26k/xmlData/data_u_t*.xml.gz"
#dataPath_p = "results/ex_03/ex_03_Re-1000.0_CFL-0.1vertices-26k/xmlData/data_p_t*.xml.gz"
#
## Load the file lists
#dataFileList_u = sorted(glob.glob(dataPath_u)) # velocity
#dataFileList_p = sorted(glob.glob(dataPath_p)) # Pressure

#===========================================================================


#===========================================================================
# Calculate lift and drag


## Pre-Allocate memory
#pressure_Drag = np.zeros(len(dataFileList_u))
##pressure_Lift = np.zeros(len(dataFileList_u))
###friction_Drag = np.zeros(len(dataFileList_u))
###friction_Lift = np.zeros(len(dataFileList_u))
##
#pressure_DragCoeff = np.zeros(len(dataFileList_u))
##pressure_LiftCoeff = np.zeros(len(dataFileList_u))
##friction_DragCoeff = np.zeros(len(dataFileList_u))
##friction_LiftCoeff = np.zeros(len(dataFileList_u))
#
#sigma_friction_Drag = np.zeros(len(dataFileList_u))
#sigma_friction_DragCoeff = np.zeros(len(dataFileList_u))
#sigma_pressure_Drag = np.zeros(len(dataFileList_u))
#sigma_pressure_DragCoeff = np.zeros(len(dataFileList_u))
#sigma_total_Drag = np.zeros(len(dataFileList_u))
#sigma_total_DragCoeff = np.zeros(len(dataFileList_u))
#
##Drag = np.zeros(len(dataFileList_u))
##DragCoeff = np.zeros(len(dataFileList_u))
## Iterate
#for iter in range(0,len(dataFileList_u)):
#    
#    print "%d of %d" % (iter+1, len(dataFileList_u))
#    
#    # Load Data
#    p = dolfin.Function(Q, dataFileList_p[iter])
#    u = dolfin.Function(V, dataFileList_u[iter])
#
#    # PRESSURE FORCES
#    
#    
#    # Pressure Lift/Drag Expression
#    #pressure_dragExpression = dolfin.inner(-p,n[0])*ds(2) # drag
#    #pressure_liftExpression = dolfin.inner(p,n[1])*ds(2) # lift
#
#    # Assembled pressure lift/drag
#    pressure_Drag[iter] = dolfin.assemble(dolfin.dot(-p,n[0])*ds(2))
#    #pressure_Lift[iter] = dolfin.assemble(dolfin.dot,(p,n[1])*ds(2))
#    
#    pressure_DragCoeff[iter] = (-pressure_Drag[iter]*2.0)/((U**2)*D)
#    #pressure_LiftCoeff[iter] = (pressure_Lift[iter]*2.0)/((U**2)*D)
#    
#    # FRICTION FORCES
#    sigma_friction_expression = 2.0*nu*epsilon(u) #- p*dolfin.Identity(u.cell().d)
#    sigma_pressure_expression = - p*dolfin.Identity(u.cell().d)
#    simga_total_expression    = 2.0*nu*epsilon(u) - p*dolfin.Identity(u.cell().d)
#    
#    #Friction Lift/Drag expression
#    sigma_friction_Drag[iter] = dolfin.assemble(dolfin.dot(dolfin.dot(sigma_friction_expression,n),eX)*ds(2))
#    sigma_pressure_Drag[iter] = dolfin.assemble(dolfin.dot(dolfin.dot(sigma_pressure_expression,n),eX)*ds(2))
#    sigma_total_Drag[iter]    = dolfin.assemble(dolfin.dot(dolfin.dot(simga_total_expression,n),eX)*ds(2))
#    #friction_dragExpression = dolfin.dot((dolfin.dot(sigma_friction, n)),eX)*dolfin.ds
#    #friction_liftExpression = (dolfin.dot(sigma_friction, n)*n[1])*dolfin.ds
#    
#    # Coefficients
#    sigma_friction_DragCoeff[iter] = (-sigma_friction_Drag[iter]*2.0)/((U**2)*D)
#    sigma_pressure_DragCoeff[iter] = (-sigma_pressure_Drag[iter]*2.0)/((U**2)*D)
#    sigma_total_DragCoeff[iter] = (-sigma_total_Drag[iter]*2.0)/((U**2)*D)
#    #friction_DragCoeff[iter] = (friction_Drag[iter]*2.0)/((U**2)*D)
#    #    friction_LiftCoeff[iter] = (friction_Lift[iter]*2.0)/((U**2)*D)
#
#        
#print 'Done.'
#===========================================================================

#===========================================================================
# Plot Lift/Drag


#py.ion()
#
#t = np.arange(0,len(dataFileList_u)*5)*((dt*U)/R)
#
#py.figure(1)
##py.plot(t[::5], friction_DragCoeff[:200],'b.-', label='friction')
#py.plot(t[::5], pressure_DragCoeff,'b-', label='pressure')
#py.axis([0,t[-1],0,2])
#py.xlabel('T')
#py.ylabel('Drag Coefficient')
#py.title('From simple pressure expression')
#py.grid()
#py.legend()
#
#
#py.figure(2)
##py.plot(t[::5], friction_LiftCoeff,'g.-', label='friction')
#py.plot(t[::5], sigma_friction_DragCoeff,'b-', label='friction')
#py.plot(t[::5], sigma_pressure_DragCoeff,'g-', label='pressure')
#py.plot(t[::5], sigma_total_DragCoeff,'k-', label='total')
#py.axis([0,t[-1],0,2])
#py.grid()
#py.xlabel('T')
#py.ylabel('Drag Coefficient')
#py.title('From $\sigma$')
#py.legend()



#py.plot(t[::5], (np.array(Drag)*2.0)/((U**2)*D),'b.-')
#py.plot(t[::5], (np.array(Lift)*2.0)/((U**2)*D),'g.-')
#===========================================================================







#===========================================================================
#===========================================================================
#===========================================================================
#===========================================================================
#===========================================================================
#===========================================================================
#===========================================================================
#===========================================================================
#===========================================================================
#===========================================================================
#===========================================================================
#===========================================================================


 #    u = dolfin.Function(V, data)
    
    #    vorticity = dolfin.project(dolfin.curl(u), Q)
    #    grad_vorticity = dolfin.project(dolfin.nabla_grad(vorticity), X)
    #    
    #    F_pressure = (-R*nu)*dolfin.inner(grad_vorticity,n)*dolfin.ds
    #    F_friction = (R*nu)*vorticity*dolfin.ds
    #    
    #    temp_F_Pressure = dolfin.assemble(F_pressure, mesh=mesh, exterior_facet_domains=bodyRegion)
    #    temp_F_Friction = dolfin.assemble(F_friction, mesh=mesh, exterior_facet_domains=bodyRegion)
    #    
    #    F_Friction.append(temp_F_Friction), F_Pressure.append(temp_F_Pressure)

#py.plot(t[::5], Cd,'.-')
#py.axis([0, t[-1], 0, 2.0])
#py.xlabel('T')
#py.ylabel('C_d')
#py.title('Evolution of Drag Coefficient')
#    #Cd = (Drag*2.0/( (U**2)*(np.pi*((D/2.0)**2))))
#    #Cl = (Lift*2.0/( (U**2)*(np.pi*((D/2.0)**2))))
#py.plot(t[::5], (np.array(F_Friction)*2.0)/((U**2)*D),'b.-')
#py.plot(t[::5], (np.array(F_Pressure)*2.0)/((U**2)*D),'g.-')

# Lift/Drag
#drag = -p*n[0]*dolfin.ds # drag
#lift = p*n[1]*dolfin.ds # lift

#vorticity = dolfin.project(dolfin.curl(u), Q)
#grad_vorticity = dolfin.project(dolfin.nabla_grad(vorticity), X)
#
#F_pressure = (-R*nu)*dolfin.dot(grad_vorticity,n)*dolfin.ds
#F_friction = (R*nu)*vorticity*dolfin.ds

# Calculate Lift/Drag function
#def calcForces():
#    Drag = dolfin.assemble(drag, mesh=mesh, exterior_facet_domains=bodyRegion)
#    #Lift = dolfin.assemble(lift, mesh=mesh, exterior_facet_domains=bodyRegion)
#    
#    Cd = (Drag*2.0/((U**2)*D))
#    Cl = (Lift*2.0/((U**2)*D))
#
#    return Lift, Drag, Cl, Cd

#def calcCoefficients(Lift, Drag):
#    
#    # Lift/Drag
#    #Lift, Drag = calcForces()
#    
#    # Lift Coefficient, Drag Coefficent [STANDARD]
#    #Cd = (Drag*2.0/( (U**2)*(np.pi*((D/2.0)**2))))
#    #Cl = (Lift*2.0/( (U**2)*(np.pi*((D/2.0)**2))))
#    
#    # Lift Coefficient, Drag Coefficent [Koumoutsakos]
#    Cd = (Drag*2.0/((U**2)*np.pi*(D/2)**2))
#    Cl = (Lift*2.0/((U**2)*np.pi*(D/2)**2))
#
#    return Cl, Cd

#===========================================================================
# Define the geometry and fluid domains

## Note: The following variables carry the location of the fluid mesh and 
## boundary function files. The fluid mesh is generated using 'Gmsh' and was
## converted to 'xml.gz' using dolfin-convert. The 'mesh' is simply the mesh
## of the fluid domain around the cylinder. The 'boundaryDomains' is a mesh
## function file that carry the identification of the domain boundaries. This
## should have been pre-defined when meshing in 'Gmsh'. The boundary [2] is the
## no-slip boundary and the boundary [3] is the outer boundary of the fluid
## domain where the far-field boundary conditions are prescribed.
#
## Define the location where simulation mesh file is saved
#mesh = "geometry/cylinder2D_Re-100_vertices-14k_MPI.xml.gz"     # mesh of the simulation fluid
#                                             # around the geometry of interest.
#                                                            
## Define the location where domain boundary function file is saved.
#boundaryDomains = "geometry/cylinder2D_Re-100_vertices-14k_MPI_facet_region.xml.gz"    # mesh boundary I.Ds
#===========================================================================


#===========================================================================
# Initialize the Navier-Stokes solver

## Initialize the NSSolver. The input parameters are the mesh 'mesh', boundary mesh
## function 'boundaryDomains', viscosity 'nu', maximum velocity 'Umax' and 
## the solver scheme type 'chorin'. The NSSolver.__init__ initializes the 
## mesh functions that are needed for inputing the initial velocity, pressure
## conditions.
#nssolver = NSSolver(mesh, boundaryDomains, nu, Umax, CFL, solverType)
#
## Initial conditions: intial velocity u0 and initial pressure p0
#u0 = dolfin.interpolate(dolfin.Constant((0.0,0.0)), nssolver.V) # ux, uy
#p0 = dolfin.interpolate(dolfin.Constant(0.0), nssolver.Q) # p
#
## Inputing the initial conditions into the solver. This will assign the u0
## and the p0 inside the solver.
#nssolver.InitialCondition(u0, p0)
#===========================================================================


#===========================================================================
# Plotting the variables

# Plot the initial conditions. NSSolver has a 'Plot' function where it is able
# to plot any plottable parameters such as u0, u1, p0, p1 and the mesh. Plots
# the variables stored in nssolver with the proper argument name.
#nssolver.Plot('mesh','u0','p0', interactive=True) # input any number of arguments with the proper variables
#===========================================================================



#===========================================================================
# Save parameters

#parameterFilePath = 'results/ex_03/info_ex_03_Re-' + str(Re) + '_CFL-' + str(CFL) + 'vertices-14k.txt'
#
#saveStyle = '%15s%30s%30s\n'
#with open(parameterFilePath,'wr+') as parameterFile:
#    parameterFile.write(saveStyle %('Variable', 'Description', 'Data'))
#    parameterFile.write(80*'-'+'\n')
#    parameterFile.write(saveStyle %('Re',   'Reynolds Number', Re))
#    parameterFile.write(saveStyle %('D',    'Charactersitic length', D))    
#    parameterFile.write(saveStyle %('U',    'Charactersitic Velocity', U))
#    parameterFile.write(saveStyle %('nu',   'Kinetmatic Viscosity', nu))
#    parameterFile.write(saveStyle %('Umax', 'Maximum Velocity', Umax))
#    parameterFile.write(saveStyle %('T',    'Time Length', T))
#    parameterFile.write(saveStyle %('CFL',  'Stability parameter', CFL))
#    parameterFile.write(saveStyle %('solverType',  'Navier-Stokes solver', solverType))
#    parameterFile.write(saveStyle %('dt',  'Time step size', nssolver.dt))
#    parameterFile.write(saveStyle %('h', 'Minimum cell size', nssolver._h))
#    parameterFile.write(saveStyle %('num_vertices', 'Number of mesh vertices', nssolver.mesh.num_vertices()))
#    parameterFile.write(saveStyle %('mesh', 'Mesh FilePath\t\t', mesh))
#    parameterFile.write(saveStyle %('boundaryDomains', 'Mesh boundary FilePath\t\t', boundaryDomains))
#
#parameterFile.close()
#===========================================================================



#===========================================================================
# Single time stepping

## Determine the time step size. For now, we wish to used the nssolver's maximum
## time step size 'dtMax'
#dt = nssolver.dtMax
#
## Feeding the velocity solution into the solver. The solver will extract
## the dirichlet velocity b.c where the solver's mesh boundary is.    
#nssolver._Step_test(U)
#
## Plot the solutions of the new time step, velocity, pressure and vorticity
#nssolver.Plot('vorticity', interactive=True)
#===========================================================================


#===========================================================================
# Storing the solutions

## .xml files for post-processing and post-calculations
## Store solution at every time-step. This save solution can be later used
## for further calculations. The 'Save_XMLDATA' saves the velocity and 
## pressure solution of the new time step u1, p1 in '.xml.gz' compressed
## format to the any location.
#nssolver.Save_XMLData(iter*nssolver.dtMax,'./results/xmlData/')
#
## .vtu/.pvd files for plotting
## Save the Plot data for further plotting in advanced plotting application
## such as 'PARAVIEW'. For such evaluations only plot of every 50-100
## iterations is needed.
#nssolver.Save_Plot('./results/')
#===========================================================================


#===========================================================================
# Iterative time stepping with NSSolver

## Determine the time step size. For now, we wish to used the nssolver's maximum
## time step size 'dtMax'
#dt = nssolver.dtMax
#
## Determine the maximum number of time steps. Number of iteration is
## determined by the number of files
#iterMax = 10# int(tend/dt)
#
#Lift = np.zeros(iterMax)
#Drag = np.zeros(iterMax)
#
#py.ion()
#
#import time
## Iterating through all the time steps. From 1 to maximum iterations.
#for iter in range(0, iterMax):
#        
#    
#    # Solving for every step    
#    nssolver._Step_test(U)
#
#    # Plot the solutions of the new time step, velocity, pressure and vorticity
#    #nssolver.Plot('vorticity') # Note: interactive=True only after the simulation run
#        
#    T_current = (iter*dt + dt)*U/(D*0.5)
#    # Display the iteratio parameters: current time, the nth iteration, run percentage.
#    print "t = %f. \t %d of %d. \t %f%% done." % (T_current, iter+1, iterMax, 100*(float(iter+1)/iterMax))
#
#    # Store solution at every time-step. This save solution can be later used
#    # for further calculations. The 'Save_XMLDATA' saves the velocity and 
#    # pressure solution of the new time step u1, p1 in '.xml.gz' compressed
#    # format
#    #if iter % 50 == 0:
#    #nssolver.Save_XMLData(T_current, 'results/ex_03/ex_03_Re-' + str(Re) + '_CFL-' + str(CFL) + 'vertices-14k/xmlData/')
#    
#    # Save the Plot data for further plotting in advanced plotting application
#    # such as 'PARAVIEW'. For such evaluations only plot of every 50-100
#    # iterations is needed.
#    #if iter % 2 == 0:
#    #    nssolver.Save_Plot(T_current, 'results/ex_03/ex_03_Re-' + str(Re) + '_CFL-' + str(CFL) + 'vertices-14k/')
#    
#    #OMEGA[iter] = 0.5*dolfin.assemble(dolfin.inner(dolfin.nabla_grad(nssolver._vort), dolfin.nabla_grad(nssolver._vort))*dolfin.dx)
#        
#    
#    # Calculate lift/drag
#    Lift[iter], Drag[iter] = calcForces() 
#    
#    # Plot data
#    py.figure(1)
#    py.scatter(iter, Lift[iter])
#        
#    py.draw()
#    time.sleep(0.05)
#    py.figure(2)
#    py.scatter(iter, Drag[iter])
#    
#    py.draw()
#    time.sleep(0.05)
#    py.show()
#    
#    
#    
#
##dolfin.interactive()
#===========================================================================


"""
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
"""
