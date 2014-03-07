#-*- coding: utf-8 -*-
__doc__ = """Navier-Stokes solver: Incremental Pressure-Correction Scheme.

Description
-----------

The solve the linear system of equations DOLFIN (PETSc) Krylov solver is
used with absolute and relative tolerance of 1E-25 and 1E-12 respectively.
To solve the velocity system, the equations where solved using GMRES
with ILU preconditioning [1]_ [2]_. 

The projection step (Darcy problem)  for :math:`u_h^n, u_h^n` is as follows,

.. math::

    \frac{u_h^n-u_h^*}{k_n} + \nabla p_h^n = 0,
    
    \nabla \cdot \left(u_h^n\right) = 0
    
which reduces to the poisson problem,
   
.. math::

    -\Delta p_h^n = \frac{-\nabla \cdot \left( u_h^*\right)}{k_n}                                                     
    
References
----------
.. [1] Valen-sendstad, B. K., Logg, A., Mardal, K., Narayanan, H., & Mortensen, 
        M. (2012). Automated Solution of Differential Equations by the Finite 
        Element Method. (A. Logg, K.-A. Mardal, & G. Wells, Eds.), 
        84. doi:10.1007/978-3-642-23099-8
            Chapter 21: A comparison of finite element schemes for the 
            incompressible Navier–Stokes equations.
                Section 21.3.1: Chorin's Projection Method
    
.. [2] Program structure:
       Kristian Valen-Sendstad <kvs@simula.no, 
       Anders Logg <logg@simula.no> 

    
:First Added:   2013-12-18
:Last Modified: 2014-02-21
:Copyright:     Copyright (C) 2014 Lento Manickathan, **pHyFlow**
:Licence:       GNU GPL version 3 or any later version 
"""

# External module
import dolfin

# Import solver base class
from pHyFlow.eulerian.base import solverBase

__all__ = ['ipcs']

class ipcs(solverBase):
    r"""
    Navier-Stokes IPCS solver class for solving the navier-stokes incompressible
    laminar problem.
        
    Usage
    -----
    .. code-block:: python
    
        NSSolver = ipcs(mesh,boundaryDomains,nu,cfl,uMax)
                 
    Parameters
    ----------
    mesh : str
           the mesh data filename (.xml.gz) generated by a FE mesh generateor 
           (GMSH) and converted to XML format.
           
    boundaryDomains : str
                      the facet boundary domain mesh data file location. The 
                      facet boundary file should be marked (int format) 
                      according to .. py:module:: pHyFlow.navierStokes.nsOptions
                      as shown below:
                                
                      (1) : Fluid Domain
                      (2) : No-Slip boundary
                      (3) : External dirichlet boundary
                      (4) : [optional] Pressure outlet.
    
                      * Note: Pressure outlet is optional.
                      
    nu : float
         the fluid kinematic viscosity :math:`\nu`.
         
    cfl : float
          the :math:`CFL` stability parameter. For explicit time-stepping,
          :math:`CFL\le1`  
    
    uMax : float
           the maximum fluid velocity :math:`U_{max}`.
    
    Attributes
    ----------
    a1 : ufl.form.Form
         the LHS of the tentative velocity expression
    
    a2 : ufl.form.Form
         the LHS of the pressure correction expression
         
    a3 : ufl.form.Form
         the LHS of the velocity correction expression
    
    A1 : dolfin.cpp.la.Matrix
         the assembled LHS of the tentative velocity expression
         
    A2 : dolfin.cpp.la.Matrix
         the assembled LHS of the pressure correction expression
         
    A3 : dolfin.cpp.la.Matrix
         the assembled LHS of the velocity correction expression.
         
    aVort : ufl.form.Form
            LHS of the variational form of vorticity            
            
    AVort : dolfin.cpp.la.Matrix
            Assembled LHS of the variation form of vorticity
            
    bcExt : dolfin.fem.bcs.DirichletBC
            The dirichlet velocity b.c at the external boundary domain (3).
    
    bcPressure : list
                 The list of pressure boundary condition.
                 * Note: empty is there is no pressure outlet.
                 
    bcNoSlip : dolfin.fem.bcs.DirichletBC
               The dirichlet velocity b.c at the no-slip boundary wall (2).
               
    bcVelocity : list
                 List of dirichlet boundary conditions.
                 
    beta : dolfin.functions.constant.Constant
           term to remove boundary stress term if problem is periodic.
                                         
    boundary_DOFCoordinates : numpy.ndarray(float64), shape (2,nDOFs)
                              array of external dirichlet boundary coordinates.
    
    boundary_VectorDOFIndex : numpy.ndarray(int64), shape (2,nDOFs)
                              array of indices of the external dirichlet 
                              boundary DOF in the velocity vector function 
                              space.
    
    boundaryDomains : dolfin.cpp.mesh.MeshFunctionSizet
                      mesh function defining the boundary domains in the mesh.
                      
    bVort : ufl.form.Form
            RHS of the variational form of vorticity.
            
    cmLocal : dolfin.cpp.mesh.Point
              the local reference point of rotation.
              
    deltaT : float
             the time step size.
             
    deltaTMax : float
                the maximum allowable time step size.
                
    f : dolfin.functions.constant.Constant
        the RHS of the navierstokes problem (source terms)
    
    F1 : ufl.form.Form
         the tentative velocity expression.    
       
    hmin : float
           the minimum mesh cell size. 
    
    k : dolfin.functions.constant.Constant
        the time step size dolfin constant.
    
    L1 : ufl.form.Form
         the RHS of the tentative velocity expression.
         
    L2 : ufl.form.Form
         the RHS of the pressure correction expression.
         
    L3 : ufl.form.Form
         the RHS of the velocity correction expression.
         
    mesh : dolfin.cpp.mesh.mesh
           the fluid mesh class
    
    n : ufl.geometry.FacetNormal
        the geometry facet normal
    
    p : dolfin.functions.function.Argument
        the navierstokes variational from pressure trial function
    
    p0 : dolfin.functions.function.Function
         the old pressure field
    
    p1 : dolfin.functions.function.Function
         the current pressure field
         
    q : dolfin.functions.function.Argument
        the pressure term test function.
        
    Q : dolfin.functions.functionspace.FunctionSpace
        the pressure function space
    
    solverName : str
                 the name of the solver. 
                 
    u : dolfin.functions.function.Argument
        the ns variational form velocity trial function
        
    U : ufl.tensors.ComponentTensor
        the linear interpolation of the old and the current velocity.         
        
    u0 : dolfin.functions.functions.Function
         the old velocity field
         
    u1 : dolfin.functions.functions.Function
         the current velocity field
         
    u1Boundary : dolfin.functions.function.Function
                 the current dirichlet boundary velocity field
                 
    uMax : float
           the maximum fluid velocity
           
    v : dolfin.functions.functions.Argument
        the velocity test function
        
    V : dolfin.functions.functionspace.VectorFunctionSpace
        the velocity vector function space
    
    w : dolfin.functions.function.Function
        the vorticity trial function
        
    x : dolfin.functions.functions.Argument
        the vorticity test function
        
    X : dolfin.functions.functionspace.FunctionSpace
        the vorticity function space
    
        
    :First Added:   2013-12-18
    :Last Modified: 2014-02-21
    :Copyright:     Copyright (C) 2014 Lento Manickathan, **pHyFlow**
    :Licence:       GNU GPL version 3 or any later version        
    """
    def __init__(self,mesh,boundaryDomains,nu,cfl,uMax):

        
        # Call the parent init (the init in solverBase)
        super(ipcs,self).__init__(mesh,boundaryDomains,nu,cfl,uMax)
        
        self.solverName = 'ipcs'
        
        # Define function space
        #self.DG = dolfin.FunctionSpace(self.mesh, "DG", 0)

        # Remove boundary stress term if problem is periodic.
        self.beta = dolfin.Constant(1) # Problem for us is not periodic
        
        # Tentative velocity step
        self.U  = 0.5*(self.u0 + self.u)
        self.F1 = (1/self.k)*dolfin.inner(self.v, self.u - self.u0)*dolfin.dx \
                    + dolfin.inner(self.v, dolfin.grad(self.u0)*self.u0)*dolfin.dx \
                    + dolfin.inner(self.epsilon(self.v), self.sigma(self.U, self.p0, self.nu))*dolfin.dx \
                    + dolfin.inner(self.v, self.p0*self.n)*dolfin.ds - self.beta*self.nu*dolfin.inner(dolfin.grad(self.U).T*self.n, self.v)*dolfin.ds \
                    - dolfin.inner(self.v, self.f)*dolfin.dx
        self.a1 = dolfin.lhs(self.F1)
        self.L1 = dolfin.rhs(self.F1)

        # Pressure correction
        self.a2 = dolfin.inner(dolfin.grad(self.q), dolfin.grad(self.p))*dolfin.dx
        self.L2 = dolfin.inner(dolfin.grad(self.q), dolfin.grad(self.p0))*dolfin.dx - (1/self.k)*self.q*dolfin.div(self.u1)*dolfin.dx

        # Velocity correction
        self.a3 = dolfin.inner(self.v, self.u)*dolfin.dx
        self.L3 = dolfin.inner(self.v, self.u1)*dolfin.dx - self.k*dolfin.inner(self.v, dolfin.grad(self.p1 - self.p0))*dolfin.dx

        # Assemble matrices
        self.A1 = dolfin.assemble(self.a1)
        self.A2 = dolfin.assemble(self.a2)
        self.A3 = dolfin.assemble(self.a3)
        
        

    def solve(self):
        """
        Solve the problem to the next time step using the dirichlet
        boundary condition.
        
        Usage
        -----
        .. code-block:: python
        
            solve()
            
        Parameters
        ----------
        None
                       
        Returns
        -------
        None returned.
        
        Attributes
        ----------
        p0 : dolfin.functions.function.Function
             the old pressure field
    
        p1 : dolfin.functions.function.Function
             the current pressure field
        
        u0 : dolfin.functions.functions.Function
             the old velocity field
         
        u1 : dolfin.functions.functions.Function
             the current velocity field
                                          
        :First Added:   2013-07-26        
        :Last Modified: 2014-02-21
        :Copyright:     Copyright (C) 2014 Lento Manickathan, **pHyFlow**
        :Licence:       GNU GPL version 3 or any later version    
        """
        
        # Previously computed solution is now the old pressure, velocity fields u0 <- u1, p0 <- p1
        self.u0.assign(self.u1) # replace u0 with u1  
        self.p0.assign(self.p1) # replace p0 with p1 

        # Compute tentative velocity step
        b = dolfin.assemble(self.L1) # assemble the RHS
        [bc.apply(self.A1, b) for bc in self.bcVelocity] # apply the velocity b.c
        dolfin.solve(self.A1, self.u1.vector(), b, "gmres", "default") # solve (u1)

        # Pressure correction
        b = dolfin.assemble(self.L2)
        if len(self.bcPressure) == 0: dolfin.normalize(b)
        [bc.apply(self.A2, b) for bc in self.bcPressure]
        dolfin.solve(self.A2, self.p1.vector(), b)
        if len(self.bcPressure) == 0: dolfin.normalize(self.p1.vector())

        # Velocity correction
        b = dolfin.assemble(self.L3)
        [bc.apply(self.A3, b) for bc in self.bcVelocity]
        dolfin.solve(self.A3, self.u1.vector(), b, "gmres", 'default')
        
        # Now that the velocity field changed we have to recalculate the 
        # vorticity
        self.recalculateVorticityFlag = True
        