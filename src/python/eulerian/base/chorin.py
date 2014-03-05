#-*- coding: utf-8 -*-
__doc__ = """ Navier-Stokes solver: Chorin's projection method.

DESCRIPTION:
    Chorin's project method (a.k.a non-incremental pressure correction
    scheme) by Chorin and Temam. The methodology is the compute the tentative
    velocity by first neglecting the pressure in the momentum equation. Then
    project the velocity on the space of divergence free vector fields. Note:
    As the velocity correct step is implemeted as the solution of a linear
    systems, discrete incompressibility constrain is not satisfied exactly. 
    But, the dirichlet b.c for velocity is applied strongly due to velocity
    correction step.
    
    The solve the linear system of equations DOLFIN (PETSc) Krylov solver is
    used with absolute and relative tolerance of 1E-25 and 1E-12 respectively.
    To solve the velocity system, the equations where solved using GMRES
    with ILU preconditioning. 
    
        Projection step (Darcy problem for u_h^n, u_h^n):
    
             ^n    ^*               
            u_h - u_h     __  ^n
            ---------  +  \/ p_h  =  0,
               k_n 
                          
                      __     ^n
                      \/ . (u_h)  =  0
                                       
        reduces to the poisson problem:
       
                                      __     ^*
                          .   ^n    - \/ . (u_h)
                       - /_\ p_h  =   ----------
                                          k_n
                                                         
    
ALGORITHM:
    1) Compute a tentative velocity (us), with proper boundary conditions for
        velocity.
    2) Compute the corrected pressure (p1), with proper boundary conditions 
        for pressure.
    3) Compute the corrected velocity (u1) using the corrected pressure (p1),
        with proper boundary condition for velocity.

REFERENCE:
    Valen-sendstad, B. K., Logg, A., Mardal, K., Narayanan, H., & Mortensen, 
    M. (2012). Automated Solution of Differential Equations by the Finite 
    Element Method. (A. Logg, K.-A. Mardal, & G. Wells, Eds.), 
    84. doi:10.1007/978-3-642-23099-8
        Chapter 21: A comparison of finite element schemes for the 
        incompressible Navier–Stokes equations.
            Section 21.3.1: Chorin's Projection Method
    
    Program structure:
        Anders Logg <logg@simula.no>, 
        Kent-Andre Mardal <kent-and@simula.no>
        
        and 
        Supervisor: Artur Palha <a.palha@tudelft.nl>
    
            
Author      :   %s
First added :   %s           
Copyright   :   %s
Licence     :   %s        
"""

__author__      = "Lento Manickathan <l.manickathan@student.tudelft.nl>"
__date__        = "2013-07-15"
__copyright__   = "Copyright (C) 2013 " + __author__
__license__     = "GNU GPL version 3 or any later version"
__doc__         %= (__author__, __date__, __copyright__ , __license__)


# Additional modules
import dolfin
from pHyFlow.eulerian.base import solverBase

__all__ = ['chorin']

class chorin(solverBase):
    def __init__(self,mesh,boundaryDomains,nu,cfl,uMax):
        """      
        ### WORK IN PROGRESS ####
        
            Initialize the solver specific parameters (for the chorin scheme).
            
            Usage
            -----
                Solver.initialize(mesh, nu, dt, f)
                
            Parameters
            ----------
                mesh    :: the geometry mesh data file location.
                ----       (type: dolfin.cpp.Mesh; single object)
                
                nu      :: the fluid kinematic viscosity.
                --         (type: dolfin.functions.constant; singe object)
                
                dt      :: the time step size.
                --         (type: dolfin.functions.constant; singe object)
                
                f       :: the Right-Hand side of the Navier-Stokes equation,
                -          i.e the source terms. Default = (0.0,0.0).
                           (type: dolfin.functions.constant; singe object)
                           
            Returns
            -------
                V       :: the vector function space of the mesh representing
                -          a vector-valued finite element function space.
                           (type: dolfin.functions.functionspace.VectorFunctionSpace; single object)
                
                Q       :: the function space of the mesh representing a 
                -          finite element function space.
                           (type: dolfin.functions.functionspace.FunctionSpace; single object)
                           
                u0      :: the class representing the initial/old velocity field
                --         in the vector function space V. Default = zeros.
                           (type: dolfin.functions.function.Function; single object)
                           
                u1      :: the class representing the new computed velocity field
                --         function in the vector function space V. Default = zeros.
                           (type: dolfin.functions.function.Function; single object)
                           
                p0      :: the class representing the intial/old pressur field
                --         in the function space Q. Default = zeros.
                           (type: dolfin.functions.function.Function; single object)
                
                p1      :: the class representing the new computed pressure field
                --         in the function space Q. Default = zeros.
                           (type: dolfin.functions.function.Function; single object)
                        
            First added: 2013-07-15
        """
        
        # Call the parent init
        super(chorin,self).__init__(mesh,boundaryDomains,nu,cfl,uMax)
              
        self.solverName = 'chorin'
        
        # Define function space
        #       The function spaces that are used for solving are the 
        #       vector function space 'V' (with 2 degree) for the velocity and
        #       function space 'Q' (with 1 degree) for the pressure. Both the
        #       function spaces belongs to the 'Lagrange'/'Continuous Galerkin' 
        #       element family.
        #        self.V = dolfin.VectorFunctionSpace(mesh, "CG", 2)
        #        self.Q = dolfin.FunctionSpace(mesh, "CG", 1)

        # Define test and trial functions
        #       Initializing the test and the trial functions for the velocity,
        #       pressure respectively.
        #        self.v = dolfin.TestFunction(self.V)
        #        self.q = dolfin.TestFunction(self.Q)
        #        self.u = dolfin.TrialFunction(self.V)
        #        self.p = dolfin.TrialFunction(self.Q)
        
        # Define functions
        self.us = dolfin.Function(self.V)   # Tentative velocity 
        #        self.u0 = dolfin.Function(self.V)   # Initial velocity before the time step
        #        self.u1 = dolfin.Function(self.V)   # Calculated velocity after the time step
        #        self.p0 = dolfin.Function(self.Q)   # Initial pressure before the time step
        #        self.p1 = dolfin.Function(self.Q)   # Calculated pressure after the time step         
        #        self.k  = dt                        # time-step
        #        self.f  = f                         # RHS (source-terms) [Default = 0]
        #        self.nu = nu                        # fluid viscosity
        
        # Tentative velocity step
        self.F1 = (1/self.k)*dolfin.inner(self.v, self.u - self.u0)*dolfin.dx \
                    + dolfin.inner(self.v, dolfin.grad(self.u0)*self.u0)*dolfin.dx \
                    + self.nu*dolfin.inner(dolfin.grad(self.v), dolfin.grad(self.u))*dolfin.dx \
                    - dolfin.inner(self.v, self.f)*dolfin.dx
        self.a1 = dolfin.lhs(self.F1)
        self.L1 = dolfin.rhs(self.F1)
        
        # Poisson problem for the pressure
        self.a2 = dolfin.inner(dolfin.grad(self.q), dolfin.grad(self.p))*dolfin.dx
        self.L2 = -(1/self.k)*self.q*dolfin.div(self.us)*dolfin.dx
        
        # Velocity update
        self.a3 = dolfin.inner(self.v, self.u)*dolfin.dx
        self.L3 = dolfin.inner(self.v, self.us)*dolfin.dx - self.k*dolfin.inner(self.v, dolfin.grad(self.p1))*dolfin.dx
        
        # Assemble matrices
        self.A1 = dolfin.assemble(self.a1)
        self.A2 = dolfin.assemble(self.a2)
        self.A3 = dolfin.assemble(self.a3)
        
        #return self.V, self.Q, self.u0, self.u1, self.p0, self.p1
        
     
    def solve(self, bcu, bcp):
        """
            Solve the problem to the next time step using the dirichlet
            boundary condition. (Note: bcp is empty).
            
            Usage
            -----
                Solver.solve(bcu, bcp)
                
            Parameters
            ----------
                bcu     :: the no-slip boundary condition on the body boundary (2)
                ---        and the dirichlet boundary condition on mesh exterior
                           boundary (3).
                           (type: )
                           
                bcp     :: the empty pressure boundary condition. Default = [].
                ---        (type:)
                           
            Returns
            -------
                u0      :: the class representing the initial/old velocity field
                --         in the vector function space V.
                           (type: dolfin.functions.function.Function; single object)
                           
                u1      :: the class representing the new computed velocity field
                --         function in the vector function space V.
                           (type: dolfin.functions.function.Function; single object)
                           
                p0      :: the class representing the intial/old pressur field
                --         in the function space Q.
                           (type: dolfin.functions.function.Function; single object)
                
                p1      :: the class representing the new computed pressure field
                --         in the function space Q.
                           (type: dolfin.functions.function.Function; single object)
                           
            First added: 2013-07-15                           
        
        """

        # Previously computed solution is now the old pressure, velocity fields u0 <- u1, p0 <- p1
        self.u0.assign(self.u1) # replace u0 with u1  
        self.p0.assign(self.p1) # replace p0 with p1  
        
        # Compute tentative velocity (us)
        b = dolfin.assemble(self.L1)                                # assemble the RHS
        [bc.apply(self.A1, b) for bc in bcu]                        # apply the velocity b.c
        #dolfin.solve(self.A1, self.us.vector(), b, 'gmres', 'ilu')  # solve (us)
        #dolfin.solve(self.A1, self.us.vector(), b, 'gmres', 'hypre_euclid')  # solve (us)
        dolfin.solve(self.A1, self.us.vector(), b, 'gmres', 'default')  # solve (us)

        
        # Compute corrected pressure (p1)
        b = dolfin.assemble(self.L2)                                # assemble the RHS
        if len(bcp) == 0: dolfin.normalize(b)                       # bcp is always 0, only dirichlet velocity b.c
        [bc.apply(self.A2, b) for bc in bcp]                        # apply the pressure b.c
        dolfin.solve(self.A2, self.p1.vector(), b)                  # solve (p1)
        #dolfin.solve(self.A2, self.p1.vector(), b, 'cg', 'hypre_amg') #TODO: remove??
        if len(bcp) == 0: dolfin.normalize(self.p1.vector())        # bcp is alway 0. only dirichlet velocity b.c    
        
        # Compute corrected velocity (u1)
        b = dolfin.assemble(self.L3)                                # assemble the RHS
        [bc.apply(self.A3, b) for bc in bcu]                        # apply the velocity b.c
        #dolfin.solve(self.A3, self.u1.vector(), b, 'gmres', 'ilu')  # solve (u1)
        dolfin.solve(self.A3, self.u1.vector(), b, 'gmres', 'default')  # solve (u1)


               
#    def modify_dt(self):
#        """
#        If :math:`\\Delta t` is modified, need to pre-assemble the LHS of
#        *tentative velocity* expression.
#        """
#        start = time.time()
#        self.A1 = dolfin.assemble(self.a1)
#        print "Time to assemble A1: " + str(time.time()-start)
#        print 'a1 assembled to A1'