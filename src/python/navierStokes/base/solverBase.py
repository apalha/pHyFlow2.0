"""
Navier-Stokes solver base class
"""

# External modules
import dolfin       # FEniCS/DOLFIN

# Import packages
import boundary

# Solver parameters
dolfin.parameters["form_compiler"]["cpp_optimize"]          = True
dolfin.parameters["krylov_solver"]["absolute_tolerance"]    = 1e-25
dolfin.parameters["krylov_solver"]["relative_tolerance"]    = 1e-12
dolfin.parameters["krylov_solver"]["monitor_convergence"]   = False
dolfin.parameters["allow_extrapolation"]                    = True
dolfin.set_log_active(False)

__all__ = ['solverBase','fluidID','noSlipID','extID']

# Navier-Stokes solver constants
fluidID  = 1
noSlipID = 2
extID    = 3

class solverBase(object):
    r"""
    Navier-Stokes solver base class
        
    Usage
    -----
    ...
                                    
                 
    Parameters
    ----------
    ...       

    Assigns
    -------
    ...    
                                  
    Returns
    -------
    ...
              
    
    :First Added:   2013-12-18
    :Last Modified: 2013-12-19
    :Copyright:     Copyright (C) 2013 Lento Manickathan, Artur Palha, **pHyFlow**
    :Licence:       GNU GPL version 3 or any later version        
    """
    
    def __init__(self,mesh,boundaryDomains,nu,cfl,uMax):
        """
        ## WORK IN PROGRESS ##

        Initialize the problem.
        """
        
        # Get the mesh from the mesh location
        self.mesh = dolfin.Mesh(mesh)
            
        # Get the boundaryDomain from the location
        self.boundaryDomains =  dolfin.MeshFunction('size_t', self.mesh, boundaryDomains)

        # Determine boundary coordinates        
        self.boundary_DOFCoordinates, self.boundary_VectorDOFIndex = boundary.locate_boundaryDOFs(self.mesh,self.boundaryDomains,3)     
        
        # Solver global parameters
        self.hmin = dolfin.MPI.min(self.mesh.hmin()) # Minimum mesh size
    
        # Calculate the maximum dt (that satisfies the CFL condition)
        # :math:`\Delta t_{max} = CFL \cdot \frac{{\Delta h_{min}}^2}{U_{max} \cdot \left( \nu + \Delta h_{min} \cdot U_{max}\right)}`
        self.dtMax = (cfl*self.hmin**2) / (uMax*(nu + self.hmin*uMax))            
        
        # Define function spaces
        self.V = dolfin.VectorFunctionSpace(self.mesh,'CG',2)# Velocity: vector function space
        self.Q = dolfin.FunctionSpace(self.mesh,'CG',1)      # Pressure: scalar function space 
        self.X = dolfin.FunctionSpace(self.mesh,'CG',1)      # Vorticity: scalar function space
        
        # Define test and trial functions
        self.u = dolfin.TrialFunction(self.V)
        self.v = dolfin.TestFunction(self.V)

        self.p = dolfin.TrialFunction(self.Q)
        self.q = dolfin.TestFunction(self.Q)

        self.w = dolfin.TrialFunction(self.X)        
        self.x = dolfin.TestFunction(self.X)
        
        # Define functions
        self.u0 = dolfin.Function(self.V)
        self.u1 = dolfin.Function(self.V)
        
        self.p0 = dolfin.Function(self.Q)
        self.p1 = dolfin.Function(self.Q)
        
        self.k  = dolfin.Constant(self.dtMax)
        self.f  = dolfin.Constant((0.,0.)) # Right-Hand side (Source terms)
        self.nu = dolfin.Constant(nu)
    
        # Define the variational problem for vorticity
        self.aVort = dolfin.inner(self.w,self.x)*dolfin.dx # LHS
        self.bVort = dolfin.inner(dolfin.curl(self.u1),self.q)*dolfin.dx # RHS
                 
        # Assemble the matrix
        self.AVort = dolfin.assemble(self.aVort) # vorticity
                
        # Define the vorticity function
        self.w = dolfin.Function(self.X)
        
        # Define boundary conditions
        self.u1Boundary = dolfin.Function(self.V)
        self.bcNoSlip = dolfin.DirichletBC(self.V, dolfin.Constant((0.,0.)), self.boundaryDomains, noSlipID)
        
        # Define the pressure boundary conditions
        self.bcPressure = []
        
        
    
    def modify_dt(self):
        r"""
        Update the time-step when modified. If :math:`\Delta t` is modified,
        the LHS of the *tentative velocity* expression needs to be updated.
        
        Usage
        -----
        ...
        
        Parameters
        ----------
        dt
        
        Assigns
        -------
        ...

        :First Added:   2013-12-18
        :Last Modified: 2013-12-18
        :Copyright:     Copyright (C) 2013 Lento Manickathan, Artur Palha, **pHyFlow**
        :Licence:       GNU GPL version 3 or any later version                
        """
        # Re-assemble the LHS
        self.A1 = dolfin.assemble(self.a1)
        
        
    def solve(self):
        """
        Solve the problem.
        """
        pass
    
    
    def vorticity(self, reCalculate=True):
        r"""
        
        ## WORKING PROGRESS ##
        
        
        Calculate the vorticity of the current calculated velocity :math:`u_1`,
        using pre-assembled matrices. The problem was solved using ``GMRES``
        with ``ILU`` preconditioning.

        Vorticity :math:`\\omega` is defined as:
                               
            :math:`\\omega = \\nabla \\times u,`
        
        and so, the variational problem becomes:
                       
            :math:`\\int\\limits_V \\, \\omega \\cdot v \\, \mathrm{d}x = \\int\\limits_V \\, \\nabla \\times u \\cdot v \\, \mathrm{d}x`
                       
        where a_vort is the LHS, b_vort is the RHS. To speed up the 
        calculation, the LHS can be preassembled (a_vort -> A_vort) as this 
        does not change and so only the RHS needs to be assembled.
               
        In order to pre-assemble the vorticity variational equations, the set
        the input argument: ``'preAssemble' = True``. This will just
        pre-assemble all the matricies for the calculation of the vorticity.
        This will be already done during the ``__init__`` of the 
        :py:class::`NavierStokes`, therefore, this process can be ignored.
        
        Usage
        -----
        .. code:block:: python
        
            vorticity()
            vorticity(preAssemble=True)*
                    
                    
            * only done in NavierStokes.__init()
               
        Parameter
        ---------
        *preAssemble* : bool
                        the boolean status to determine if vorticity problem
                        should be pre-assembled
        
        
        Assigns
        -------
        W : dolfin.functions.functionspace.FunctionSpace
            the finite element function space to project the curl of the 
            velocity (vorticity)

        _vort : dolfin.functions.function.Function
                the test function of the variational problem to the vorticity.
                           
        _A_vort : dolfin.cpp.la.Matrix
                  the assembled LHS of the variational problem, PETScMatrix. 
                
        _b_vort : ufl.form.Form
                  the RHS of the variational problem.
               
                           
        Returns
        -------
        _vort : dolfin.functions.function.Function
                the solution to the variation problem, i.e the vorticity of 
                the current calculated velocity.
                           
                           
        :First Added:   2013-07-15
        :Last Modified: 2013-12-10
        :Copyright:     Copyright (C) 2013 Lento Manickathan, Artur Palha, **pHyFlow**
        :Licence:       GNU GPL version 3 or any later version    
                   
        """
        
        #        # Pre-Assemble the vorticity matricies
        #        if preAssemble is True:
        #
        #                # Function Space for vorticity
        #                self.W = dolfin.FunctionSpace(self.mesh,'CG',1)
        #            
        #                # Trial functions            
        #                self.w = dolfin.TrialFunction(self.W)#(self.Q) # vorticity trial function
        #                self.q = dolfin.TestFunction(self.W)#(self.Q)
        #                
        #                # Variational problem for vorticity
        #                self.a_vort = dolfin.inner(self.w,self.q)*dolfin.dx #LHS
        #                self.b_vort = dolfin.inner(dolfin.curl(self.u1),self.q)*dolfin.dx #RHS
        #                 
        #                # Assemble the matrix
        #                self.A_vort = dolfin.assemble(self.a_vort) # vorticity
        #                
        #                # Define the vorticity function
        #                self.w = dolfin.Function(self.W)
        #                
        #        else:
        # Calculate vorticity if it has not been calculated
        if reCalculate:
            # Compute Vorticity
            b = dolfin.assemble(self.bVort) # Assemble RHS
            dolfin.solve(self.AVort, self.w.vector(), b, 'gmres', 'default') # solve for vorticity
        
        # Else: Vorticity has already been calculated [just return vorticity]
        return self.w        



    def boundaryConditions(self,vx,vy):
        r"""
        Extract
        """
        
        # Transfer the vertex values of the functions
        self.u1Boundary.vector()[self.boundary_VectorDOFIndex[0]] = vx
        self.u1Boundary.vector()[self.boundary_VectorDOFIndex[1]] = vy
        
        # Exterior boundary condition
        self.bcExt = dolfin.DirichletBC(self.V, self.u1Boundary, self.boundaryDomains, extID)
        
        # Collect all the boundary conditions
        self.bcVelocity = [self.bcNoSlip, self.bcExt]
