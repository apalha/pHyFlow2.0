#-*- coding: utf-8 -*-
__doc__ = """
Navier-Stokes solver base class
"""

__all__ = ['solverBase']

# External modules
import dolfin       # FEniCS/DOLFIN
import numpy

# Import pHyFlow packages
#from pHyFlow.aux.customDecorators import simpleGetProperty

# Import eulerian functions
from pHyFlow.eulerian.base import boundary
from pHyFlow.eulerian import eulerianOptions

# Solver parameters
dolfin.parameters["form_compiler"]["cpp_optimize"]          = eulerianOptions.FORM_COMPILER['cpp_optimize']
dolfin.parameters["krylov_solver"]["absolute_tolerance"]    = eulerianOptions.KRYLOV_SOLVER['absolute_tolerance']
dolfin.parameters["krylov_solver"]["relative_tolerance"]    = eulerianOptions.KRYLOV_SOLVER['relative_tolerance']
dolfin.parameters["krylov_solver"]["monitor_convergence"]   = eulerianOptions.KRYLOV_SOLVER['monitor_convergence']
dolfin.parameters["allow_extrapolation"]                    = eulerianOptions.ALLOW_EXTRAPOLATION
dolfin.set_log_active(eulerianOptions.SET_LOG_ACTIVE)


class solverBase(object):
    r"""
    Navier-Stokes solver base class for solving the navier-stokes incompressible
    laminar problem.

    Usage
    -----
    .. code-block:: python

        <customSolver>(solverBase)

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

    hmin : float
           the minimum mesh cell size.

    k : dolfin.functions.constant.Constant
        the time step size dolfin constant.

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

    u : dolfin.functions.function.Argument
        the ns variational form velocity trial function

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

        # Get the mesh from the mesh location
        self.mesh = dolfin.Mesh(mesh)

        # Store the local mesh coordinates
        self.mesh_localCoordinates = self.mesh.coordinates().copy()

        # Get the boundaryDomain from the location
        self.boundaryDomains =  dolfin.MeshFunction('size_t', self.mesh, boundaryDomains)

        self.cmLocal = dolfin.Point() # (0.,0.,0.)

        # Determine boundary coordinates
        #self.boundary_DOFCoordinates, self.boundary_VectorDOFIndex = boundary.locate_boundaryDOFs(self.mesh,self.boundaryDomains,3)
        # Solver global parameters
        self.hmin = dolfin.MPI.min(self.mesh.mpi_comm(),self.mesh.hmin()) # Minimum mesh size

        # Mesh components
        self.n  = dolfin.FacetNormal(self.mesh) # normal function
        self.eX = dolfin.Constant((1.0, 0.0)) # cartesian unit x-vector
        self.eY = dolfin.Constant((0.0, 1.0)) # cartesian unit y-vector
        self.ds = dolfin.Measure("ds")[self.boundaryDomains] # boundary line integrator

        # Calculate the maximum dt (that satisfies the CFL condition)
        # :math:`\Delta t_{max} = CFL \cdot \frac{{\Delta h_{min}}^2}{U_{max} \cdot \left( \nu + \Delta h_{min} \cdot U_{max}\right)}`
        self.cfl = cfl
        self.uMax = uMax
        self.deltaTMax = (self.cfl*self.hmin**2) / (self.uMax*(nu + self.hmin*uMax))

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

        self.k  = dolfin.Constant(self.deltaTMax)
        self.f  = dolfin.Constant((0.,0.)) # Right-Hand side (Source terms)
        self.nu = dolfin.Constant(nu)

        # Define the variational problem for vorticity
        self.aVort = dolfin.inner(self.w,self.x)*dolfin.dx # LHS
        #self.bVort = dolfin.inner(dolfin.curl(self.u1),self.q)*dolfin.dx # RHS
        self.bVort = dolfin.inner(dolfin.curl(self.u1),self.x)*dolfin.dx # RHS

        # Assemble the matrix
        self.AVort = dolfin.assemble(self.aVort) # vorticity

        # Define the vorticity function
        self.w = dolfin.Function(self.X)

        # Define boundary conditions
        self.u1Boundary = dolfin.Function(self.V)
        self.bcNoSlip = dolfin.DirichletBC(self.V, dolfin.Constant((0.,0.)), self.boundaryDomains, eulerianOptions.ID_NOSLIP_BOUNDARY)

        # Determine boundary indices and boundary coordinates
        #self.boundary_VectorDOFIndex = boundary.vectorDOF_boundaryIndex(self.V,self.boundaryDomains,options.ID_EXTERNAL_BOUNDARY)
        self.boundary_VectorDOFIndex = boundary.vectorDOF_boundaryIndex(self.mesh,self.boundaryDomains,eulerianOptions.ID_EXTERNAL_BOUNDARY,2)
        #self.boundary_DOFCoordinates = boundary.vectorDOF_coordinates(self.mesh,self.V,self.boundary_VectorDOFIndex[0])

        # Define the pressure boundary conditions
        if eulerianOptions.ID_PRESSURE_OUTLET in self.boundaryDomains.array():
            # If the boundary domain has pressure outlet
            self.bcPressure = [dolfin.DirichletBC(self.Q, dolfin.Constant(0.), self.boundaryDomains, eulerianOptions.ID_PRESSURE_OUTLET)]
        else:
            # If there is no pressure outlet
            self.bcPressure = []


    def set_deltaT(self,deltaT):
        r"""
        Update the time-step when modified. If :math:`\Delta t` is modified,
        the LHS of the *tentative velocity* expression needs to be updated.

        Usage
        -----
        .. code-block:: python

            set_deltaT(deltaTNew)

        Parameters
        ----------
        deltaT : float
                 the new time step size. The deltaT should satisfy
                 :math:`\Delta t \le \Delta t_{max}`

        Returns
        -------
        None returned.

        Attributes
        ----------
        deltaT : float
                 the new time step size. The deltaT should satisfy
                 :math:`\Delta t \le \Delta t_{max}`

        k : dolfin.functions.constant.Constant
            the time step size dolfin constant.

        A1 :

        :First Added:   2013-12-18
        :Last Modified: 2014-02-21
        :Copyright:     Copyright (C) 2013 Lento Manickathan, Artur Palha, **pHyFlow**
        :Licence:       GNU GPL version 3 or any later version
        """
        # Check if deltaT should be set to max
        if deltaT == 'max':
            self.deltaT = self.deltaTMax

        # check is delta smaller than deltaTMax
        elif deltaT > self.deltaTMax:
            raise ValueError("""Time step 'deltaT = %g' is larger than the
                             maximum allowable time step 'deltaTMax = %g', at
                             CFL = %g""" % (deltaT, self.deltaTMax, self.cfl))

        else:
            self.deltaT = deltaT

        # Assign the delta T to k
        self.k.assign(self.deltaT)

        # Re-assemble the LHS
        self.A1 = dolfin.assemble(self.a1)


    def solve(self):
        "Solve the problem ."
        pass


    def vorticity(self):
        r"""

        Calculate the vorticity of the current calculated velocity :math:`u_1`,
        using pre-assembled matrices. The problem was solved using ``GMRES``
        with ``ILU`` preconditioning.

        Vorticity :math:`\omega` is defined as:

            :math:`\omega = \nabla \times u,`

        and so, the variational problem becomes:

            :math:`\int\limits_V \, \omega \cdot v \, \mathrm{d}x = \int\limits_V \, \nabla \times u \cdot v \, \mathrm{d}x`

        where a_vort is the LHS, b_vort is the RHS. To speed up the
        calculation, the LHS can be preassembled (a_vort -> A_vort) as this
        does not change and so only the RHS needs to be assembled.

        Usage
        -----
        .. code:block:: python

            vorticity()

        Parameter
        ---------
        recalculate : bool
                      should the vorticity be recalculated

        Returns
        -------
        w : dolfin.functions.function.Function
            the solution to the variation problem, i.e the vorticity of
            the current calculated velocity.

        Attributes
        ----------
        w : dolfin.functions.functionspace.FunctionSpace
            the solution to the variation problem, i.e the vorticity of
            the current calculated velocity.

        :First Added:   2013-07-15
        :Last Modified: 2014-02-21
        :Copyright:     Copyright (C) 2014 Lento Manickathan, **pHyFlow**
        :Licence:       GNU GPL version 3 or any later version

        """

        # Calculate vorticity if it has not been calculated
        if self.recalculateVorticityFlag:
            # Compute Vorticity
            b = dolfin.assemble(self.bVort) # Assemble RHS
            dolfin.solve(self.AVort, self.w.vector(), b, 'gmres', 'default') # solve for vorticity

        # Else: don't recalculate vorticity

        # Return vorticity field
        return self.w



    def boundaryConditions(self,vx,vy):
        r"""
        Function to apply the dirichlet velocity boundary condition at the
        external boundary domain definied by
        .. py:module:: pHyFlow.navierStokes.nsOptions.ID_EXTERNAL_BOUNDARY.

        Usage
        -----
        .. code-block:: python

            boundaryConditions(vx,vy)

        Parameters
        ----------
        vx : numpy.ndarray(float64), shape (nDOFs,)
             the :math:`x` component of the dirichlet velocity boundary
             condition at the navierstokes DOF boundaries.

        vy : numpy.ndarray(float64), shape (nDOFs,)
             the :math:`y` component of the dirichlet velocity boundary
             condition at the navierstokes DOF boundaries.

        Returns
        -------
        None returned.

        Attributes
        ----------
        bcExt : dolfin.fem.bcs.DirichletBC
                The dirichlet velocity b.c at the external boundary domain (3).

        bcPressure : list
                     The list of pressure boundary condition.
                     * Note: empty is there is no pressure outlet.

        bcNoSlip : dolfin.fem.bcs.DirichletBC
                   The dirichlet velocity b.c at the no-slip boundary wall (2).

        bcVelocity : list
                     List of dirichlet boundary conditions.

        u1Boundary : dolfin.functions.function.Function
                     the current dirichlet boundary velocity field

        :First Added:   2013-07-15
        :Last Modified: 2014-02-21
        :Copyright:     Copyright (C) 2014 Lento Manickathan, **pHyFlow**
        :Licence:       GNU GPL version 3 or any later version
        """

        # Transfer the vertex values of the functions
        self.u1Boundary.vector()[self.boundary_VectorDOFIndex[0]] = vx
        self.u1Boundary.vector()[self.boundary_VectorDOFIndex[1]] = vy

        # Exterior boundary condition
        self.bcExt = dolfin.DirichletBC(self.V, self.u1Boundary, self.boundaryDomains, eulerianOptions.ID_EXTERNAL_BOUNDARY)

        # Collect all the boundary conditions
        self.bcVelocity = [self.bcNoSlip, self.bcExt]


    def updateMeshPosition(self, cmGlobal, thetaLocal):
        r"""
        Function to update the position of the mesh.

        Usage
        -----
        .. code-block:: python

            updateMeshPosition(cmGlobal, thetaLocal)

        Parameters
        ----------
        cmGlobal : numpy.ndarray(float64), shape (2,)
                   the :math:`x,y` position of the mesh local reference
                   point (0.,0.) in the global coordinates.

        thetaLocal : float, unit (rad)
                     the local rotational angle :math:`\theta` of the mesh
                     domain. Therefore, the rotation will be done about local
                     reference point (0.,0.), i.e cmGlobal in the global
                     coordinate system.

        Returns
        -------
        None

        Attributes
        ----------
        mesh : dolfin.cpp.mesh.mesh
               the fluid mesh class at the new global position

        :First Added:   2014-03-06
        :Last Modified: 2014-03-06
        :Copyright:     Copyright (C) 2014 Lento Manickathan, **pHyFlow**
        :Licence:       GNU GPL version 3 or any later version

        """
        # Reset the mesh coordinates (global) to local coordinates
        self.mesh.coordinates()[:] = self.mesh_localCoordinates

        # Rotate the mesh by thetaLocal around (0.,0.)
        self.mesh.rotate(numpy.rad2deg(thetaLocal),2,self.cmLocal) # cmLocal = [0.,0.]

        # Move the mesh to new global position cmGlobal.
        self.mesh.coordinates()[:] += cmGlobal #


    def epsilon(self, u):
        "Return symmetric gradient."
        return 0.5*(dolfin.grad(u) + dolfin.grad(u).T)


    def sigma(self, u, p, nu):
        "Return stress tensor."
        return 2.0*nu*self.epsilon(u) - p*dolfin.Identity(u.cell().d)


    def FrictionalForces(self):
        r"""
        Calculates the friction force components acting on the
        body (no-slip boundary).
        """

        # Friction force : x-component
        Fx = - (dolfin.assemble(dolfin.dot(dolfin.dot((2.0*self.nu*self.epsilon(self.u1)),
                                                      self.n),
                                           self.eX)*self.ds(eulerianOptions.ID_NOSLIP_BOUNDARY)))

        # Friction force : y-component
        Fy = - (dolfin.assemble(dolfin.dot(dolfin.dot((2.0*self.nu*self.epsilon(self.u1)),
                                                      self.n),
                                           self.eY)*self.ds(eulerianOptions.ID_NOSLIP_BOUNDARY)))

        # Return the friction forces
        return numpy.array([Fx, Fy])


    def PressureForces(self):
        r"""
        Calculate the pressure-gradient force components acting on the
        body (no-slip boundary).
        """

        # Pressure force: x-component
        Fx = - (dolfin.assemble(dolfin.dot(dolfin.dot(-self.p1*dolfin.Identity(self.u1.cell().d),
                                                      self.n),
                                           self.eX)*self.ds(eulerianOptions.ID_NOSLIP_BOUNDARY)))

        # Pressure force: x-component
        Fy = - (dolfin.assemble(dolfin.dot(dolfin.dot(-self.p1*dolfin.Identity(self.u1.cell().d),
                                                      self.n),
                                           self.eY)*self.ds(eulerianOptions.ID_NOSLIP_BOUNDARY)))

        # return the pressure forces
        return numpy.array([Fx, Fy])


    def Forces(self):
        r"""
        Calculates the total forces acting on the body
        """

        # Total force: x-component (drag)
        Fx = -(dolfin.assemble(dolfin.dot(dolfin.dot((2.0*self.nu*self.epsilon(self.u1) - self.p1*dolfin.Identity(self.u1.cell().d)),
                                                     self.n),
                                          self.eX)*self.ds(eulerianOptions.ID_NOSLIP_BOUNDARY)))

        # Total force: y-compnent (Lift)
        Fy = -(dolfin.assemble(dolfin.dot(dolfin.dot((2.0*self.nu*self.epsilon(self.u1) - self.p1*dolfin.Identity(self.u1.cell().d)),
                                                     self.n),
                                          self.eY)*self.ds(eulerianOptions.ID_NOSLIP_BOUNDARY)))

        # return the total forces
        return numpy.array([Fx, Fy])


    def totalCirculation(self):
        """
        Function to calculate the total circulation
        """
        return dolfin.assemble(self.vorticity()*dolfin.dx)

    def vectorDOF_boundaryCoordinates(self):
        r"""
        vectorDOF_boundaryCoordinates : numpy.ndarray(float64), shape (2,nVectorBoundaryDOFs)
                                        the :math:`x,y` boundary coordinates of the
                                        vector DOFs.
        """
        # Extract the boundary coordinates from complete vector dof coordiantes
        return numpy.copy(boundary.vectorDOF_coordinates(self.mesh,self.V)[:,self.boundary_VectorDOFIndex[0]])


    def vectorDOF_coordinates(self):
        r"""
        vectorDOF_coordinates : numpy.ndarray(float64), shape (2,nVectorDOFs)
                                the :math:`x,y` coordinates of the vector DOFs.
        """
        # return the copy of the reshaped, unrepeated (unique) vector DOFs coordinates
        return numpy.copy(self.V.dofmap().tabulate_all_coordinates(self.mesh).reshape(-1,2).T[:,::2])
