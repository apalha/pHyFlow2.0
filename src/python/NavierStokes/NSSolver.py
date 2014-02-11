#-*- coding: utf-8 -*-
__doc__ = """ Navier-Stokes solver for time-stepping with Dirichlet Boundary condition.
DESCRIPTION:
    Navier-Stokes Finite Element Method Solver for incompressible fluid. 
    Also, a turbulence scheme is ot implemented and therefore the problem
    should be laminar. Solves the problem for a single time-step with 
    appropriate initial velocity, pressure field and final velocity dirichlet 
    boundary condition. The solver is based on the FEniCS/DOLFIN 1.2.0 
    Navier-Stokes solvers and is implemented based on the algorithms provided
    by their example cases.
    
REFERENCE:
    Solver:
        FEniCS/DOLFIN 1.2.0
            http://fenicsproject.org/about/            
        
    Book:    
        Valen-sendstad, B. K., Logg, A., Mardal, K., Narayanan, H., & Mortensen, 
        M. (2012). Automated Solution of Differential Equations by the Finite 
        Element Method. (A. Logg, K.-A. Mardal, & G. Wells, Eds.), 
        84. doi:10.1007/978-3-642-23099-8
            Chapter 21: A comparison of finite element schemes for the 
            incompressible Navier–Stokes equations.    
   
   Algorithm:
       NSbench is a collection of Python scripts based on FEniCS/DOLFIN for 
       benchmarking a number of common methods for solving the incompressible 
       Navier-Stokes equations. Anders Logg <logg@simula.no> - NSBENCH
           https://launchpad.net/nsbench
           
           Implemented solver algorithms:
               
               - Chorin projection scheme ['chorin']

Author      :   %s
First added :   %s           
Copyright   :   %s
Licence     :   %s  
      
"""

__author__      = "Lento Manickathan <l.manickathan@student.tudelft.nl> & Artur Palha <a.palha@tudelft.nl>"
__date__        = "2013-07-15"
__copyright__   = "Copyright (C) 2013 " + __author__
__license__     = "GNU GPL version 3 or any later version"
__doc__         %= (__author__, __date__, __copyright__ , __license__)

__all__ = ["NSSolver"]

# Additional modules
import dolfin       # FEniCS/DOLFIN
import numpy        # Standard numerical package
import fenicstools  # Fenicstools - postprocessing tools for FEniCS by mikaem
from solvers import Solver # Solver schemes
import os           # Library of OS routines

# Solver parameters
dolfin.parameters["form_compiler"]["cpp_optimize"]          = True
dolfin.parameters["krylov_solver"]["absolute_tolerance"]    = 1e-25
dolfin.parameters["krylov_solver"]["relative_tolerance"]    = 1e-12
dolfin.parameters["krylov_solver"]["monitor_convergence"]   = False
dolfin.parameters["allow_extrapolation"]                    = True
dolfin.set_log_active(False)
#dolfin.parameters["num_threads"] = 6 # Multi-threading 



class NSSolver:
    __doc__ = __doc__
    """
        Reviews:
        -------
            (1) [2013-07-13], with Artur Palha:
                - Class restructured.
                - removed problem.py, moved all the functions to here.
                
            (2) [2013-07-17], with Artur Palha:
                - Implement the velocity array input for dirichlet b.c at
                  exterior domain instead of the velocity function.
                  
            (3) [2013-08-16]
                - Structured interpolation of vorticity added.
                - '_interpolateVorticity_initialize' : Initializes the function
                - 'interpolateVoriticy': probes the vorticity at the dof of
                    the interpolation grid
                    
            (4) [2013-08-20]
                - aribitrary number of substeps possible now.
                - Implemented in the .step function.

    """
    
    def __init__(self, mesh, boundaryDomains, nu, uMax, CFL, noSlipBC=True, solverType='chorin', interpolationMesh_bounds=None, interpolationMesh_density=None):
        """
            Navier-Stokes solver using FEniCS/DOLFIN Finite Element Solver.
        
            Usage
            -----
                NSSolver = NSSolver(mesh, boundaryDomains, nu, uMax, solverType='chorin')
                
            Parameters
            ----------
                mesh    :: the geometry mesh data file location. Mesh generated
                ----       using Gmsh (finite element mesh generator) or similar
                           finite element generator, converted to DOLFIN XML
                           format (.xml, .xml.gz).
                           (type: str; shape: singe word)
                    
                boundaryDomains :: the facet boundary domain mesh data file
                ---------------    location. facet boundary domain data file
                                   generated using 'dolfin-convert', format (.xml, .xml.gz).
                                   (type: str; shape: single word)
                                  
                nu      :: the fluid kinematic viscosity.
                --         (type: float (64bits); shape: single value)
                         
                uMax    :: the maximum fluid velocity to compute the CFL number.
                ----       (type: float (64bits); shape: single value)
                
                CFL     :: the simulation CFL number. For explicit time
                ---        marching, CFL should be less than 1. 
                           (type: float (64bits); shape: single value)
                           
                noSlipBC :: (optional) the boolean status whether or not if there is a 
                            no-slip boundary (i.e geometry) is found in the 
                            simulation. In other words, if boundary domain (2)
                            exists. Default is True.
                            (type: bool; shape: single value)
                           
                solverType :: (optional) the solver algorithm name to be used for
                              solving the NS problem. Default='chorin'.
                              (type: str; shape: single word)
                           
                interpolationMesh_bounds  :: (optional) the mesh bounds for the structured 
                                             interpolation grid. The vorticity 
                                             will be transfered to the blobs from 
                                             this structured grid. The meshBounds
                                             should be of form: [xmin, ymin, xmax, ymax].
                                             (type: float (64bits); shape: [4,])
                
                
                interpolationMesh_density :: (optional) the mesh density for the structured 
                                             interpolation grid. The vorticity
                                             will be tranfered to the blobs from
                                             this grid. The meshDensity should be 
                                             of form: number of vertices in [Nx, Ny].
                                             (type: int; shape: [2,])
                   
            Returns
            -------
                NSSolver :: solver with pre-initialized matrices for solving the 
                            navier-stokes problem using the chosen solver algorithm,
                            default is chorin project scheme
                            (type: instance)
                  
                  
            First added:     2013-07-15

            Copyright (C) 2013 Lento Manickathan, Artur Palha
            
        """
        
        # Parameter for storing data solution files
        self._uFile = None 
        self._pFile = None
        self._wFile = None
        
        self._saveDir = None
    
        # Load the mesh and boundary Domain mesh for user input
        self._LoadMesh(mesh, boundaryDomains)

        # Defining the boundary identifications
        self._noSlipBC          = noSlipBC
        self._noSlipDomainID    = 2 # no-slip region
        self._exteriorDomainID  = 3 # exterior domain
        
        # Solver global parameters
        self._cfl   = CFL   # CFL number
        self._nu    = nu    # Fluid viscosity
        self._uMax  = uMax  # Maximum fluid velocity
        self._h     = dolfin.MPI.min(self.mesh.hmin())  # Minimum mesh size
    
        # Calculate the maximum dt (that satisfies the CFL condition)
        self._calc_dtMax()
    
        # Setting up Navier-Stokes constants in dolfin format
        self.f  = dolfin.Constant((0.0,0.0))    # Right-Hand side (Source terms)
        self.nu = dolfin.Constant(self._nu)     # Fluid viscosity (in DOLFIN format)
        self.dt = dolfin.Constant(self.dtMax)   # Time step
    
        # Choose the solver: (Chorin or other solvers)
        # - Various Navier-stokes solving algorithms can be imported. All the
        # algorithms should be stored in './solvers'. The solver is then
        # imported in as part of this module into -> 'self.solver'.
        self.Solver = Solver(solverType)
            
        # Initialize the solver parameters
        # Pre-allocate all the variables for the FEM problem
        self.V, self.Q, self.u0, self.u1, self.p0, self.p1 = self.Solver.initialize(self.mesh, self.nu, self.dt, self.f)

        # Pre-assemble matricies for calculating vorticity
        self._vorticityCalculated = False   # Vorticity not calculated yet
        self.vorticity(preAssemble=True)    # Assemble the vorticity
        
        # Initialize the InterpolateVorticity
        if interpolationMesh_bounds and interpolationMesh_density is not None:
            self._interpolation_meshBounds = interpolationMesh_bounds   # Define the bounds of the interpolation mesh
            self._interpolation_meshDensity = interpolationMesh_density # Defin the mesh resolution
            self._interpolateVorticity_initialize()                     # Initialize the interpolation of vorticity
        
        # Find the boundary notes and store them
        self._boundaryDOFs_evaluate() # find the boundary DOFs
        self._u0_boundary = dolfin.Function(self.V) # Initializing u1_boundary function 
        self._u1_boundary = dolfin.Function(self.V) # Initializing u0_boundary function
        self._bcu_grad    = self._u0_boundary.vector().array() # Initializing bcu gradient array

        # Initializing the no-slip bc.
        if self._noSlipBC:
            self._g2 = dolfin.Constant((0.0,0.0))
            self._bc2 = dolfin.DirichletBC(self.V, self._g2, self.boundaryDomains, self._noSlipDomainID) 
        
        # Pressure b.c is None
        self.bcp = [] # Note: problem does not have pressure b.c
        
        
    def boundaryConditions(self, u1_boundary):
        '''
            Extracts dirichlet velocity boundary condition from the target
            velocity field using appropriate boundary domain i.d and with 
            appropriate function space.
            
                Boundary I.D:
                    - 2 :: no-slip boundary
                    - 3 :: exterior mesh boundary
                    
            'u1_boundary' can be either; type(None) [Default] if exterior boundary
            condition is zero, type(DOLFIN.functions.functions.Function) or 
            type(NUMPY.ndarray). If the data-type is later, then it will be
            applied accordingly.
        
            Usage
            -----
                bcu, bcp = NSSolver.boundaryConditions(u1_boundary)
            
            Parameters
            ----------
                u1_boundary :: (optional) the target velocity field. The velocity 
                               field should have the valid dirichlet boundary
                               condition at the exterior mesh boundary.
                               (type: dolfin.functions.function.Function; single object    .or
                                type: numpy.ndarray (float64), size: size of vector DOF    .or
                                type: NoneType, size: single object )
                        
            Assigns
            -------
                bcu :: the collection of all the dirichlet velocity boundary
                ---    condition. No-slip and exterior b.c respectively.
                       (type: list(dolfin.fem.bcs.DirichletBC), size: [2]) 
                
                bcp :: the collection of all the dirichlet pressure boundary
                ---    condition. No dirichlet pressure condition for solver
                       (type: empty list, size: [])
                       
           Returns
           -------
               (-)
                       
        '''         
            
        # If the input data 'u1_boundary' is a NUMPY array in the proper order as defined by the dof map
        if isinstance(u1_boundary, numpy.ndarray):
            
            # Transfer the vertex values of the function
            self._u1_boundary.vector()[self.boundary_vectorDOF[0,:]] = u1_boundary[0,:] # x component
            self._u1_boundary.vector()[self.boundary_vectorDOF[1,:]] = u1_boundary[1,:] # y component       
            
            # Exterior boundary condition
            bc3 = dolfin.DirichletBC(self.V, self._u1_boundary, self.boundaryDomains, self._exteriorDomainID)
            
         # If the input file u1_boundary is a DOLFIN function in function space V            
        elif isinstance(u1_boundary, dolfin.functions.function.Function):
            
            # Exterior boundary condition from 'u1_boundary'
            bc3 = dolfin.DirichletBC(self.V, u1_boundary, self.boundaryDomains, self._exteriorDomainID)
            
        else:
            
            # Else if not a DOLFIN FUNCTION nor a NUMPY array
            raise NameError('Input boundary conditions dolfin function or numpy array!')
            
        # Collect the b.c.s
        if self._noSlipBC:
            self.bcu = [self._bc2, bc3]
        else:
            self.bcu = [bc3]
        
        

    def initialConditions(self, u0, p0):
        '''
            Assigns the initial velocity field u0 (u_{0}) and the initial 
            pressure field p0 (p_{0}) to the problem. u_{0} function should 
            be of the 'NSSolver' vector function space (V) and p_{0} should
            be of the 'NSSolver' function space (Q). Therefore, the input
            data should be already interpolated to the NSSolver to the function
            spaces.
        
            Usage
            -----
                NSSolver.initialConditions(u0, p0)
            
            Parameters
            ----------
                u0  :: the initial/old velocity field in the vector
                --     function space V. (interpolated onto V)
                       (type: dolfin.functions.function.Function; single object)
                        
                           
                p0  :: the intial/old pressur field in the function
                --      space Q. (interpolated onto Q)
                       (type: dolfin.functions.function.Function; single object)
                       
                       
            Assigns
            -------
                u0  :: the initial/old velocity field in the vector
                --     function space V. (interpolated onto V)
                       (type: dolfin.functions.function.Function; single object)
                       
                p0  :: the intial/old pressur field in the function
                --      space Q. (interpolated onto Q)
                       (type: dolfin.functions.function.Function; single object)
                    
            Returns
            -------
                (-)
                
        '''
        
        # Check if u0 function space is the same as the initialized function space
        # The input initial velocity field, pressure field should have 
        # been already interpolated into the solver function spaces.
        if u0.function_space() != self.V:
            raise RuntimeError, 'u0 function space is not the same as the '+\
                                'initialized function space nssolver.V! '+\
                                'Interpolate to the solver function space! Idiot!'
                                
        if p0.function_space() != self.Q:
            raise RuntimeError, 'p0 function space is not the same as the ' +\
                                'initialized function space nssolver.Q! ' +\
                                'Interpolate to the solver function space! Idiot!'
        
        # Set the initial condition 
        # Note: This automatically updates self.Solver.[u0,u1,p0,p1] as references are already made
        self.u0.assign(u0), self.u1.assign(u0) # also u1 because in nssolver, u0<-u1
        self.p0.assign(p0), self.p1.assign(p0)
        
        # Initial boundary condition
        self._u0_boundary.vector()[:] = self.u0.vector()[:]
        
        
    def interpolateVorticity(self):#TODO: update doc
        """
            Interpolate vorticity onto the blobs. The vorticity is probed at the
            DOFs of the interpolation grid from the Navier-Stokes domain and is
            then used to later transfer to the blobs.
            
            Usage
            -----
                NSSolver.interpolateVorticity(blobs)
                
            Parameters
            ----------
                blobs   :: ...
                -----
                
            Assigns
            -------
                :math:`\\omega` :: vorticity at the blob locations.
                
            Returns
            -------
                (-)
            
        """
        
        # Remove the old data
        if self._interpolation_probes.number_of_evaluations() != 0:
            self._interpolation_probes.clear()
            
        # Probe the data
        self._interpolation_probes(self.vorticity())
        
        # Transfer the vorticity to blobs
        # self.transfer2blobs(self._interpolation_probes.array())
        
        # Return the data
        return self._interpolation_xGrid, self._interpolation_yGrid, self._interpolation_probes.array().reshape(self._interpolation_xGrid.shape)
        
        

    def step(self, u1_boundary, dt, nSteps=1):
        """
            Navier-Stokes solver does a single time step. Solves the problem 
            for a single time-step with appropriate target dirichlet 
            velocity boundary condition and proper time step parameter. The 
            time step size should be with in the CFL stability bound.
            
            Multi-Stepping possible by changing 'nSteps' > 1. If 'nSteps' is
            greater than 1, :math:'dt_{sub-step} = \\frac{dt}{nSteps}'.
            
            Usage
            -----
                NSSolver.step(u1_boundary, dt=dtMax)
                NSSolver.step(u1_boundary, dt=dt, nSteps=4)
            
            Parameters
            ----------
                u1_boundary     :: the target velocity field. The velocity 
                -----------        field should have the valid dirichlet 
                                   boundary condition at the exterior mesh 
                                   boundary.
                                   (type: dolfin.functions.function.Function; single object)
            
                dt  :: the time step size between the initial 
                       solutions (u0, p0) and the target solution (u1,p1). 
                       The time step size should be within the stability bound 
                       (dt <= dt_{max}) definied by the CFL condition ('cfl'). 
                       (type: float (64bits); single value)
                       
                nSteps :: (optional) the number of substep per current step. This 
                          should be chosen such that the substep time step size 
                          should be  within the stability bounds of the 'cfl' number.
                          Default : nstep = 1.
                          (type: int; single value)
                   
            Assigns
            -------
                u0  :: the class representing the initial/old velocity field
                --     in the vector function space V. Default = zeros.
                       (type: dolfin.functions.function.Function; single object)
                               
                u1  :: the class representing the new computed velocity field
                --     function in the vector function space V. Default = zeros.
                       (type: dolfin.functions.function.Function; single object)
                               
                p0  :: the class representing the intial/old pressur field
                --     in the function space Q. Default = zeros.
                       (type: dolfin.functions.function.Function; single object)
                    
                p1  :: the class representing the new computed pressure field
                --     in the function space Q. Default = zeros.
                       (type: dolfin.functions.function.Function; single object)
                       
            Returns
            -------
                (-)
                       
        """
        
        # Determine the time-step
        self._dt = dt/nSteps # Fraction to number of steps
            
        # Check if dt is within stability bounds
        self._check_dt(self._dt) # dt<=dtMax
        
        
        # If the input b.c is DOLFIN function
        if isinstance(u1_boundary, dolfin.functions.function.Function):
            
            # Calculate the gradient of the b.c            
            self._bcu_grad = (u1_boundary.vector().array() - self._u0_boundary.vector().array())/nSteps
            
        elif isinstance(u1_boundary, numpy.ndarray):
            
            # Calculate the gradient of the b.c
            self._bcu_grad[self.boundary_vectorDOF[0,:]] = (u1_boundary[0,:] - self._u0_boundary.vector()[self.boundary_vectorDOF[0,:]])/nSteps
            self._bcu_grad[self.boundary_vectorDOF[1,:]] = (u1_boundary[1,:] - self._u0_boundary.vector()[self.boundary_vectorDOF[1,:]])/nSteps
            
        else:
            raise RuntimeError("Unknown format. Either DOLFIN function or NUMPY.ndarray!")


        # Multi-Stepping: from 1 to the nth-step
        for k in range(1,nSteps+1):
            
            print "substep: " + str(k)
            
            # Calculate the b.c of the sub-step
            self._u1_boundary.vector()[:] = self._u0_boundary.vector().array() + self._bcu_grad*k
            
            # Set the boundary conditions
            self.boundaryConditions(self._u1_boundary)
            
            # Solve the problem
            self.Solver.solve(self.bcu, self.bcp)
        

        # Change the boolean status of vorticity calculation
        self._vorticityCalculated = False
        
        
        # Store the b.c of the old step with the new step
        self._u0_boundary.vector()[:] = self._u1_boundary.vector()[:]
            
            
            
    def plot(self, *args, **kwargs):
        """
            Plot any DOLFIN plottable parameters using DOLFIN.plot. This
            module is able to plot multiple parameters and a given call.
            
            Usage
            -----
                NSSolver.plot("mesh")
                NSSolver.plot("mesh", interactive=True)
                NSSolver.plot("u1")
                NSSolver.plot("u1","p1")
                NSSolver.plot("u1","p1","vorticity")
                NSSolver.plot("u1","p1","vorticity", interactive=True)
            
            Parameters
            ----------
                *args   :: the NSSolver attributes names. Attributes that 
                -----      corresponds to the DOLFIN plottable parameters such 
                           as 'u0', 'u1', 'p0', 'p1', 'mesh', 'boundaryDomains'
                           and 'vorticity'. If argument is 'vorticity', the 
                           function will compute the vorticity right then (just-in-time).
                           (type: str; shape: *)
                           
                *kwargs :: additional plotting parameters such as plot
                           interactivity. To interact set 'interactive'=True.
                           (type: dict; shape: *)
            Returns
            -------
                DOLFIN.plot :: the input parameter is plotted using the DOLFIN 
                               plot function, with title = parameter name, 
                               key = parameter name.
                               (type: dolfin.common.plotting.plot)
                            
        """
        
        # Plots the parameter in the input argument
        for parameter in args:
            # Vorticity is computed only when necessary
            if parameter is 'vorticity':
                dolfin.plot(self.vorticity(), title=parameter, key=parameter) # Pre-assembled
            else:
                # Plot other parameters with its title and key
                dolfin.plot(getattr(self,parameter), title=parameter, key=parameter)        
        
        # Turn interactive on
        if 'interactive' in kwargs:
            dolfin.interactive()
        
        
        
    def saveData(self, dataType, t, saveDir):
        """
            Saves the pressure, velocity and the vorticity data of the Navier-
            stokes domain in the appropriate data format ('dataType').
        
            If .xml.gz format:  Stores the velocity 'u1' and pressure solution 'p1' 
            files in compressed form: .xml.gz format.
            
            If .pvd format: the stored solution can then
            be plotted using '.pvd' readers such as PARAVIEW.
            
                Typically stores the files to (e.g):
                        ./results/
                
            Usage
            -----
                NSSolver.saveData('pvd', t, './results/pvdData/')
                NSSolver.saveData('xml', t, './results/xmlData/')
                
            Parameters
            ----------
                dataType    :: the solution save format. Either 'pvd' or 'xml'
                --------       for saving solution in .pvd or .xml.gz (compressed)
                               format.
                               (type: str; size: single word)            
            
                t   ::  the time value of the the solution that is saved. The 
                -       time value is used to label the .xml.gz file.
                        (type: float (64bits), size: single value)
                        
                saveDir     :: the location where to save the xml data.
                -------        (type: str; singe word)
                
            Returns
            -------
                .pvd files  :: [data_u.pvd, data_p.pvd, data_w.pvd], velocity
                ----------     pressure, and vorticity respectively. The .pvd
                               files hold the library for the .vtu data files.
                               (type: data file (.pvd), size: 3*single file)
 
                .vtu files  :: [data_u*.vtu, data_p*.vtu, data_w*.vtu], velocity
                ----------     pressure and vorticity data files at each time
                               instant. For every iteration a new data file is
                               created; the file names are iterated 000000,000001,...
                               (type: data file (*.vtu), size: 3 * (number of save iterations))
               
        """
        
        # Make save Directory
        if self._saveDir is None:
            self._saveDir = saveDir            
            os.makedirs(self._saveDir) # Make leaf directories
        
        if dataType == 'pvd':
        
            # Create files for saving
            if self._uFile is None:
                self._uFile = dolfin.File(saveDir + "data_u.pvd", "compressed")   # Pressure
            if self._pFile is None:
                self._pFile = dolfin.File(saveDir + "data_p.pvd", "compressed")   # Velocity
            if self._wFile is None:
                self._wFile = dolfin.File(saveDir + "data_w.pvd", "compressed")   # Vorticity
                
            # Every other save iteration, dolfin.File automatically saves and names document
            self._uFile << (self.u1, t)  # Pressure 
            self._pFile << (self.p1, t)  # Velocity
            self._wFile << (self.vorticity(), t) # Vorticity
            
        elif dataType == 'xml':

            # Define the .xml save file path
            uFile = dolfin.File(saveDir + "data_u_" + "t%f" % t + ".xml.gz")
            pFile = dolfin.File(saveDir + "data_p_" + "t%f" % t + ".xml.gz")
            
            # Storing the velocity, pressure solution
            uFile << self.u1 # Velocity
            pFile << self.p1 # pressure
            
        else:
            # If unknown dataType
            raise NameError("Save format 'dataType' unknown. Either 'xml' format or 'pvd'!")
  
    

    def vorticity(self, **kwargs):
        """
            Calculate the vorticity of the current calculated velocity (u1),
            using pre-assembled matrices. The problem was solved using
            GMRES with ILU preconditioning.

                Vorticity (vort) is defined as:         
                               
                   .. math ::
                               
                       \\omega = \\nabla \\times u,
                               
                and so, the variational problem becomes:
                       
                   .. math ::
                   
                       \\int\\limits_V \\, \\omega \\cdot v \\, \mathrm{d}x = \\int\\limits_V \\, \\nabla \\times u \\cdot v \\, \mathrm{d}x
                       
               
               where a_vort is the LHS, b_vort is the RHS. To speed up the
               calculation, the LHS can be preassembled (a_vort -> A_vort)
               as this does not change and so only the RHS needs to be assembled.
               
               In order to pre-assemble the vorticity variational equations,
               the set the input argument: 'preAssemble' = True. This will just
               pre-assemble all the matricies for the calculation of the vorticity.
               This will be already done during the '__init__' of the 'NSSOLVER',
               therefore, this process can be ignored.
               
                       
            Usage
            -----
                    Solver.vorticity()
                    
                    ** Solver.vorticity(preAssemble=True)
                    
                    
                    ** [Will be done in NSSOLVER.__init__]
                
            Parameters
            ----------
                vort    :: the test function of the variational problem to
                           the vorticity.
                           (type: dolfin.functions.function.Function; single object)
               
                A_vort  :: the assembled LHS of the variational problem,
                           PETScMatrix. 
                           (type: dolfin.cpp.la.Matrix; size: single object)
                
                b_vort  :: the RHS of the variational problem.
                           (type: ufl.form.Form; size: single Form)
                           
            Returns
            -------
                vort    :: the solution to the variation problem, i.e 
                           the vorticity of the current calculated velocity
                           (type: dolfin.functions.function.Function; single object)
                           
        """
        
        # Pre-Assemble the vorticity matricies
        if 'preAssemble' in kwargs:
            
                # Trial functions            
                self._vort = dolfin.TrialFunction(self.Q) # vorticity trial function
                self._q = dolfin.TestFunction(self.Q)
                
                # Variational problem for vorticity
                self._a_vort = dolfin.inner(self._vort,self._q)*dolfin.dx            #LHS
                self._b_vort = dolfin.inner(dolfin.curl(self.u1),self._q)*dolfin.dx #RHS
                
                # Assemble the matrix
                self._A_vort = dolfin.assemble(self._a_vort) # vorticity
                
                # Define the vorticity function
                self._vort = dolfin.Function(self.Q)
                
        else:
            
            # Calculate vorticity if it has not been calculated
            if self._vorticityCalculated is not True:
            
                # Compute Vorticity
                b = dolfin.assemble(self._b_vort) # Assemble RHS
                dolfin.solve(self._A_vort, self._vort.vector(), b, 'gmres', 'default') # solve for vorticity
                
                # The vorticity has been calculated
                self._vorticityCalculated = True # change overhead status
        
            # Else: Vorticity has already been calculated [just return vorticity]
            return self._vort
            
        
        
    def _boundaryDOFs_evaluate(self):
        """
            Compute the mapping from vertices to degrees of freedom only at the boundary.
            This saves time, since we do not have to loop over all the cells of the domain.
            If the for loop was replaced by a simple c++ function call, it would be much faster.
        
            Parameters
            ----------
                mesh    :: the mesh where we want to compute the mapping between 
                ----       vertices and degrees of freedom.
                           (type: dolfin.cpp.Mesh; single object)
                        
                boundaryDomains     :: A MeshFunction which assigns integer 
                ---------------        values to edges of the mesh.
                                       (type: dolfin.cpp.mesh.MeshFunctionSizet; single object) 
                                    
                exteriodDomainID    :: An integer value which we use to identify 
                ----------------       the boundaryedges we want to use. The 
                                       MeshFunction should have this value 
                                       assigned to some of the edges. 
                                       (type: int; single value)
                                       
                p   :: Order of the finite element space
                -      (type: int; single value) 
            
            Assigns
            -------
                boundary_DOFCoordinates     ::  coordinates of the dofs
                -----------------------         (type: numpy.ndarray (float64); size: ([2, number of DOFs))
                                            
                boundary_vectorDOF  :: stores the mapping from vertex number to
                ------------------     degree of freedom number for vector functions.
                                       (type: numpy.ndarray (float64); size: ([2, number of DOFs))
                                       
            Returns
            -------
                (-)
                
        """

        # Order of the finite element space
        p = 2
    
        self._boundaryMesh = dolfin.BoundaryMesh(self.mesh,'exterior')
        
        # get the connectivity of the mesh
        self.mesh.init(1) # initialize edges
        mesh_topology = self.mesh.topology() # get the mesh topology
        conn12 = mesh_topology(1,2) # get the cells to which an edge belongs
            
        # get the number of edges
        n_edges = self.boundaryDomains.size()
    
        # get the indices of the boundary edges, as given by the boundaryFunction
        boundaryEdges = numpy.arange(n_edges)
        boundaryEdges = boundaryEdges[self.boundaryDomains.array() == self._exteriorDomainID]
        
        # Define the function spaces associated to the mesh
        # for now only CG of order 1 works for sure
        V = dolfin.FunctionSpace(self.mesh, "CG", p)
        Vv = dolfin.VectorFunctionSpace(self.mesh, "CG", p, dim=2)
    
        # Get the degree of freedom maps
        dm = V.dofmap()
        dms = [Vv.sub(i).dofmap() for i in range(2)] # the i=0 refers to x component
                                                     # the i=1 refers to y component
        
        # compute the total number of dofs
        n_dofs = dm.global_dimension()
        
        # Allocate memory space of the array that stores the mapping from
        # vertex number to degree of freedom number for both scalar and vector
        # valued functions (note that for vector valued functions we need to indices
        # since the x component and the y component have different indices)    
        #vert_to_dofs = numpy.zeros(n_dofs, dtype=numpy.uintp)
        self.boundary_vectorDOF = numpy.zeros([2,n_dofs], dtype=numpy.uintp)
            
        # Allocate memory space of the arrays that store the coordinates of the dofs
        self.boundary_DOFCoordinates = numpy.zeros([2,n_dofs])
        
        # get the cells of the mesh
        #mesh_cells = mesh.cells()
        
        # Since we are only interested in the boundary vertices we just loop over
        # the boundary edges as given by the MeshFunction and stored in boundaryEdges
        # each time we loop we get the numbering of the degrees of freedom of the
        # vertices of the cell that contains the boundary edge and we construct
        # the mapping from vertex number to degree of freedom number
        for edge in boundaryEdges:
            # get the cell index to which the edge belongs
            cell_ind = conn12(edge)[0] 
           
            # get the dofs of the cells
            cell_dofs = dm.cell_dofs(cell_ind)

            # add the degree of freedom numbering and the vertex coordinates of the dofs
            self.boundary_DOFCoordinates[:,cell_dofs] = dm.tabulate_coordinates(dolfin.Cell(self.mesh,cell_ind)).T
            
            # now we do the same but for vector valued quantites
            # since, in this case, our vector have two quantities (the x and y components)
            # we have to iterate over the two sub dof maps
            for i, (dms_i, dmcs_i) in enumerate(zip(dms, dms)):
                dms_i.cell_dofs(cell_ind)
                self.boundary_vectorDOF[i,cell_dofs] = numpy.array(dms_i.cell_dofs(cell_ind),dtype=numpy.uintp)
    
        # get the list of dofs at the boundary
        scalar_dof_at_boundary = numpy.array(dolfin.DirichletBC(V, dolfin.Constant(0.0), self.boundaryDomains, self._exteriorDomainID).get_boundary_values().keys())
                
        # Assigns vertex to dof mapping                
        self.boundary_DOFCoordinates = self.boundary_DOFCoordinates[:,scalar_dof_at_boundary]
        self.boundary_vectorDOF      = self.boundary_vectorDOF[:,scalar_dof_at_boundary]
        
        
    def _calc_dtMax(self):
        """
            Calculate the maximum time step size 'dtMax' for which the 
            problem is still stable, i.e within the CFL stability bound.
            
                Equation for dt_{max}: 
                                             CFL x dh_{min}^2
                    dt_{max}    =   ------------------------------------
                                    U_{max} . ( nu + dh_{min} x U_{max} )    
        
            Usage
            -----
                NSSolver._calc_dtMax(cfl, h, uMax, nu)            
            
            Parameters
            ----------
                cfl     :: the CFL value, defined by the Courant-Friedrichs-
                ---        Lewy condition. For explicit time-marching schemes,
                           typical CFL = 1.0. Lower value ensure more stable
                           solution, but returns much smaller time step size.
                           (type: float (64bits); size: single value)
                
                h       :: the small mesh grid size. The smallest mesh size
                -          determines the stability bound. h <- mesh.hmin()
                           (type: float (64bits): size: single value)
                
                uMax    :: the maximum fluid velocity. For cylinders, typical
                ----       maximum fluid velocity is U_{max} = 2.0 x U_{inf}, 
                           from potential flow solutions.
                           (type: float (64bits): size: single value)
                           
                nu      :: the fluid kinematic viscosity.
                --         (type: float (64bits): size: single value)
                
            Assigns
            -------
                dtMax   :: the maximum allowable solver time step size. If 
                -----      beyond this parameter. The problem is beyond the 
                           stability bound.
                           (type: float (64bits): size: single value)
                       
            Returns
            -------
                (-)                     
                           
        """
        
        # Calculate maximum time step size dtMax
        self.dtMax = (self._cfl*self._h**2) / (self._uMax*(self._nu + self._h*self._uMax))
        

        
    def _check_dt(self, dt):
        """
            Used to check if the input 'dt' is within the stability range. 
            This is true if (dt<=dtMax). If the dt is larger than dtMax,
            'RuntimeError' is raised, else the input 'dt' is used by the
            solver.
            
            Usage
            -----
                NSSolver._check_dt(dt,dtMax,cfl)
                
            Parameters
            ----------
                dt  :: The time step size between the initial 
                --     solutions (u0, p0) and the target solution (u1,p1). 
                       The time step size should be within the stability bound 
                       (dt <= dt_{max}) definied by the CFL condition ('cfl'). 
                       If 'dt' is not given, 'dtMax' is chosen. Default: dt=None.
                       (type: float (64bits); single value)
                       
                dtMax   :: the maximum allowable solver time step size. If 
                -----      beyond this parameter. The problem is beyond the 
                           stability bound.
                           (type: float (64bits): size: single value)              
            
                cfl     :: the CFL value, defined by the Courant-Friedrichs-
                ---        Lewy condition. For explicit time-marching schemes,
                           typical CFL = 1.0. Lower value ensure more stable
                           solution, but returns much smaller time step size.
                           (type: float (64bits); size: single value)
                           
            Assigns
            -------
                _dt :: The time step size between the initial 
                --     solutions (u0, p0) and the target solution (u1,p1).
                       Note: the returned dt is within the stability bound
                       (type: float (64bits); single value)
            
            Returns
            -------
                (-)
                       
        """
        
        # Check if dt is larger than dtMax
        if dt > self.dtMax:
            raise RuntimeError, "Time step 'dt' is beyond stability limit," +\
                                " at CFL = %f, 'dt' should be <%f" % (self._cfl, self.dtMax)
        else:
            self._dt = dt
        
        
        
    def _interpolateVorticity_initialize(self):
        """
            Initialize the InterpolateVorticity Function. InterpolateVorticity
            will be used to interpolate the vorticity to the blobs.
        
        #TODO:update
        """
                
        # Making references
        x0,y0,x1,y1 = self._interpolation_meshBounds
        nx,ny =  self._interpolation_meshDensity 
        
        # Generate interpolation mesh
        self._interpolation_xGrid, self._interpolation_yGrid = numpy.meshgrid(numpy.linspace(x0,x1,nx),
                                                                              numpy.linspace(y0,y1,ny))
                                                                              
        # Format the coordinates                                                                              
        self._interpolation_coor = numpy.append(self._interpolation_xGrid.reshape(-1,1),self._interpolation_yGrid.reshape(-1,1),axis=1)                                                                        

        # Generate dolfin interpolation mesh
        self._interpolation_mesh = dolfin.RectangleMesh(x0,y0,x1,y1,nx,ny)

        # Generate function space
        self._interpolation_functionSpace = dolfin.FunctionSpace(self._interpolation_mesh, 'CG', 1)
        
        # Structured dofmap
        self._interpolation_dofMap = self._interpolation_functionSpace.dofmap().vertex_to_dof_map(self._interpolation_mesh)
        
        # Defining the probing
        self._interpolation_probes = fenicstools.Probes(self._interpolation_coor.flatten(), self.Q)
        
        # Structured vorticity field
        #self._interpolation_vorticity = dolfin.Function(self._interpolation_functionSpace)
    
    
    
    
    def _LoadMesh(self, mesh, boundaryDomains):
        """      
            Retrieve the mesh and boundary-domains data for the input path.
            The files will be read and loaded to the meshData and boundaryDomainsData
            
            Usage
            -----
                 NSSolver._LoadMesh(mesh, boundaryDomains)
                 
            Parameters
            ----------
                mesh    :: the geometry mesh data file location. Mesh generated
                ----       using Gmsh (finite element mesh generator) or
                           similar finite element generator, converted to 
                           DOLFIN XML format (.xml, .xml.gz).
                           (type: str; shape: singe word)
                    
                boundaryDomains     :: the facet boundary domain mesh data 
                ---------------        file location. facet boundary domain 
                                       data file generated using 'dolfin-convert', 
                                       format (.xml, .xml.gz).
                                       (type: str; shape: single word)
                                   
            Assigns
            -------
                mesh    :: the loaded mesh file in the proper DOLFIN format
                ----       (type: dolfin.cpp.Mesh; single object)
                               
                boundaryDomains     :: the loaded mesh boundary file with
                ---------------        appropriate boundary I.D. (2) is 
                                       no-slip boundary condition, and 
                                       (3) is the dirichlet velocity b.c
                                       (type: dolfin.cpp.mesh.MeshFunctionSizet; single object)
                                       
            Returns
            -------
                (-)                                       
                                           
        """
        
        # Get the mesh from the mesh location
        self.mesh = dolfin.Mesh(mesh)
        
        # Get the boundaryDomain from the location
        self.boundaryDomains =  dolfin.MeshFunction("size_t", self.mesh, boundaryDomains) # load the boundaries identification



##########################################################################################
##    Used for checking if the solver is functioning accordingly. Will be removed !!    ##        
##########################################################################################            


    
    #TODO: remove when everything is working           
    def _boundaryConditions_standard(self, U, T):
        '''
            Simple free-stream boundary condition. boundaryDomain must have 4 regions.
            
            Domains:
            
                domain 2 :: no-slip boundary
                domain 3 :: in-flow boundary
                domain 4 :: out-flow boundary
                domain 5 :: (top-bottom) free-stream boundary
                
            Parameters
            ----------
                ...
                
            Assigns
            -------
                ...
                
            Returns
            -------
                (-)
                                
        '''
        
        # Invoke pertubation (rotating cylinder): Ref: Lecointe and Piquet (1984)
        if T >= 3.0 and T<3.5:
            
            # First Pertubation
            u_r     = 0.0   # Radial velocity
            u_theta = 0.15  # Tangential velocity
            
        elif T>=3.5 and T<=5.0:
            # Second pertubation
            u_r     = 0.0
            u_theta = -0.25
        else:
            # Unpertubated flow (non-rotating cylinder)
            u_r     = 0.0
            u_theta = 0.0
            
        # Define the rotational velocity field of the cylinder
        g2  = dolfin.Expression(('u_r*cos(atan2(x[1],x[0])) - u_theta*sin(atan2(x[1],x[0]))',
                                 'u_r*sin(atan2(x[1],x[0])) + u_theta*cos(atan2(x[1],x[0]))'),
                                  u_r = u_r, u_theta=u_theta)
                                  
        bc2 = dolfin.DirichletBC(self.V, g2, self.boundaryDomains, 2)
            
        # Create inflow boundary condition
        g3  = dolfin.Constant((U,0.0))
        bc3 = dolfin.DirichletBC(self.V, g3, self.boundaryDomains, 3)
        
        # Create outflow boundary condition
        g4  = dolfin.Constant(0)
        bc4 = dolfin.DirichletBC(self.Q, g4, self.boundaryDomains, 4)
        
        # Create free stream boundary condition for velocity
        g5  = dolfin.Constant((U, 0.0))
        bc5 = dolfin.DirichletBC(self.V, g5, self.boundaryDomains, 5)
            
        # Collect boundary conditions
        self.bcu = [bc2, bc3, bc5]
        self.bcp = [bc4]
        
        

    #TODO: not in the final version     
    def _step_standard(self, U=None, dt=None, T=None):
        """
            Standard Stepping (with pressure outlet)
            
            Parameters
            ----------
                ...
                
            Returns
            -------
                (-)
                
        """

        # Set the time-step 
        if dt is None:
            self._dt = self.dtMax
        else:
            # Check if dt is within stability bound
            self._check_dt(dt)
      
        # Set the time-step
        self.dt.assign(self._dt)
        
        # Set boundary conditions
        if U is not None:
            self._boundaryConditions_standard(U, T)
        
        # Solve the problem
        self.Solver.solve(self.bcu, self.bcp)
        
        # Change boolean status of vorticity calculation
        self._vorticityCalculated = False
        
        

#    def saveXMLData(self, t, saveDir):
#        """
#            Stores the velocity 'u1' and pressure solution 'p1' files in .xml 
#            (or compressed form: .xml.gz) format. This saved solution can then
#            be loaded into DOLFIN solver using 'dolfin.Function'. NOTE: it is
#            adviced to make a separate folder for the xml datas. The number of
#            file that is generated is highly and might slow the system down 
#            if you view it with GUI file manager. 
#            
#                Stores the files to the input directiory (e.g):
#                    ./results/xmlData/
#                    
#                Note: make directory 'xmlData' in 'results' dir beforehand
#            
#            Usage
#            -----
#                NSSolver.Save_XMLData(t, saveDir)
#                NSSolver.Save_XMLData(NSSolver.dtmax*iteration, saveDir)
#            
#            Parameters
#            ----------
#                t   ::  the time value of the the solution that is saved. The 
#                -       time value is used to label the .xml.gz file.
#                        (type: float (64bits), size: single value)
#                        
#                saveDir     :: the location where to save the xml data.
#                -------        (type: str; singe word) 
#                                                            
#        
#            Returns
#            -------
#                .xml.gz files   :: [data_u_t*.xml.gz, data_p_t*.xml.gz], the
#                -------------      velocity and pressure (u1, p1) respectively.
#                                   (type: data file (*.xml.gz); size: 3 * (number of iteration))            
#            
#        """
#
#        # Define the .xml save file path
#        uFile = dolfin.File(saveDir + "data_u_" + "t%f" % t + ".xml.gz")
#        pFile = dolfin.File(saveDir + "data_p_" + "t%f" % t + ".xml.gz")
#        
#        # Storing the velocity, pressure solution
#        uFile << self.u1 # Velocity
#        pFile << self.p1 # pressure        



    


#    def Initialize_InterpolationGrid(self, dl, N, origin=[0.0, 0.0], angle=0.0):
#        '''
#            Initialize the structured Interpolation Grid. This grid will be
#            used to interpolate the parameters (such as vorticity) from the
#            Navier-Stokes domain to the vortex particle domain.
#
#            Usage
#            -----
#            
#                    
#            Parameters
#            ----------
#            
#            Assigns
#            -------
#
#            Returns
#            -------
#        
#        '''
#        
#        # Set-up the interpolation grid 2D      
#        x0, y0 = origin[0] - dl[0]/2, origin[1] - dl[1]/2
#        x1, y1 = origin[0] + dl[0]/2, origin[1] + dl[1]/2
#        nx, ny = N
#        
#        # Generate structured mesh
#        self.structured_mesh = dolfin.RectangleMesh(x0,y0,x1,y1,nx,ny)
#    
#        # Rotate Mesh
#        self.structured_mesh.rotate(angle)
#        
#        # Structured FunctionSpace
#        self.structured_Q = dolfin.FunctionSpace(self.structured_mesh, 'CG', 1)
#        
#        # Structured dofmap
#        self.structured_dofmap = self.structured_Q.dofmap().vertex_to_dof_map(self.structured_mesh)
#        
#        # Defining the probing
#        self.structured_probes = fenicstools.Probes(self.structured_mesh.coordinates().flatten(), self.Q)
#        
#        # Structured vorticity field
#        self.structured_vorticity = dolfin.Function(self.structured_Q)



#    def step(self, u1_boundary=None, dt=None): #TODO:update doc
#        """
#            Navier-Stokes solver does a single time step. Solves the problem 
#            for a single time-step with appropriate target dirichlet 
#            velocity boundary condition and proper time step parameter. The 
#            time step size should be with in the CFL stability bound.
#            
#            Usage
#            -----
#                NSSolver.step(u1_boundary, dt=dt)
#                NSSolver.step(u1_boundary, dt=None)
#            
#            Parameters
#            ----------
#                u1_boundary     :: the target velocity field. The velocity 
#                -----------        field should have the valid dirichlet 
#                                   boundary condition at the exterior mesh 
#                                   boundary.
#                                   (type: dolfin.functions.function.Function; single object)
#            
#                dt  :: (optional) the time step size between the initial 
#                       solutions (u0, p0) and the target solution (u1,p1). 
#                       The time step size should be within the stability bound 
#                       (dt <= dt_{max}) definied by the CFL condition ('cfl'). 
#                       If 'dt' is not given, 'dtMax' is chosen. Default: dt=None.
#                       (type: float (64bits); single value)
#                   
#            Returns
#            -------
#                u0  :: the class representing the initial/old velocity field
#                --     in the vector function space V. Default = zeros.
#                       (type: dolfin.functions.function.Function; single object)
#                               
#                u1  :: the class representing the new computed velocity field
#                --     function in the vector function space V. Default = zeros.
#                       (type: dolfin.functions.function.Function; single object)
#                               
#                p0  :: the class representing the intial/old pressur field
#                --     in the function space Q. Default = zeros.
#                       (type: dolfin.functions.function.Function; single object)
#                    
#                p1  :: the class representing the new computed pressure field
#                --     in the function space Q. Default = zeros.
#                       (type: dolfin.functions.function.Function; single object)
#               
#        """
#
#        # Determine the Time-step size      
#        if dt is None:
#            self._dt = self.dtMax # if dt is not given, dtMax is the dt
#        else:
#            # Check if dt is within stability bound
#            self._check_dt(dt)
#            
#        # Assign the dt to the DOLFIN 'dt' parameter
#        self.dt.assign(self._dt) # Note: also assigns in the solver
#        
#        # Set the boundary conditions
#        self.boundaryConditions(u1_boundary)
#        
#        # Solve the problem
#        self.Solver.solve(self.bcu, self.bcp)
#        
#        # Change boolean status of vorticity calculation
#        self._vorticityCalculated = False



        
#        # If the input file u1_boundary is a DOLFIN function in function space V
#        if isinstance(u1_boundary, dolfin.functions.function.Function):
#            
#            # Exterior boundary condition from 'u1_boundary'
#            bc3 = dolfin.DirichletBC(self.V, u1_boundary, self.boundaryDomains, self._exteriorDomainID) 
#        
#        # If the input data 'u1_boundary' is a NUMPY array in the proper order as defined by the dof map
#        elif isinstance(u1_boundary, numpy.ndarray):
#            
#            # Transfer the vertex values of the function
#            self._u1_boundary.vector()[self.boundary_vectorDOF[0,:]] = u1_boundary[0,:] # x component
#            self._u1_boundary.vector()[self.boundary_vectorDOF[1,:]] = u1_boundary[1,:] # y component             
#            
#            # Exterior boundary condition from '_u1_boundary' function
#            bc3 = dolfin.DirichletBC(self.V, self._u1_boundary, self.boundaryDomains, self._exteriorDomainID)
#            
#            # Collect boundary conditions
#            self.bcu = [self._bc2, bc3] # no-slip, exterior b.c
#            self.bcp = []               # Note: problem does not have pressure b.c
#          
#          
#        # No 'u1_boundary' velocity. Only no-slip boundary condition is found
#        elif u1_boundary is None:
#            
#            # Collect boundary conditions
#            self.bcu = [self._bc2]
#            self.bcp = []
#            
#        else:
#            raise NameError('Input boundary conditions dolfin function or numpy array!')
#   