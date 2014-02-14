"""
Functions related Navier-Stokes boundary
"""

# External modules
import numpy
import dolfin

def vectorDOF_boundaryIndex(mesh,boundaryDomains,boundaryID,pOrder):
    r"""
    Function to locate vector DOF (Degrees of freedom) boundary index in the
    vector function space velocity field array.
    
    The returned vector DOF boundary index **vectorDOF_boundaryIndex** can be
    used to calculate the coordinates of the vector DOF boundary at any instant.
    
    Note: The index is valid even after any transformation to the fluid domain.
    
    Usage
    -----
    .. code-block:: python
    
        xIndices,yIndices = vectorDOF_boundaryIndex(V,boundaryDomains,boundaryID)
        
    or

    .. code-block:: python
    
        xyIndices = vectorDOF_boundaryIndex(V,boundaryDomains,boundaryID)
        
    Parameters
    ----------
    mesh : dolfin.cpp.mesh
           the dolfin mesh of the fluid domain.
           
    boundaryDomains : dolfin.cpp.mesh.MeshFunctionSizet
                      A MeshFunction which assigns integer values to the 
                      edges of the mesh.
                      
    boundaryID : int64
                 The boundary ID integer value found inside the mesh function
                 **boundaryDomains**. 
                 Note: - fluid              :: 1
                       - no-slip boundary   :: 2
                       - external boundary  :: 3
                       - pressure outlet    :: 4
                       
    pOrder : int
            the order of the finite element space of the vector function space.
                       
    Returns
    -------
    xyIndices : numpy.ndarray(float64), shape (2,N_boundaryDOFs)
                the :math:`x,y` indicies of the vector boundary DOFs.
                       
    :First Added:   2014-02-14
    :Last Modified: 2014-02-14
    :Copyright:     Copyright (C) 2014 Lento Manickathan, **pHyFlow**
    :License:       GNU GPL version 3 or any later version
    """
 
    # Temporary Function space, same order as vector function space
    Q = dolfin.FunctionSpace(mesh,'CG',pOrder)
    
    # Extract the boundary index values   
    iIndices = numpy.array(dolfin.DirichletBC(Q,dolfin.Constant(0.),
                                              boundaryDomains,boundaryID).get_boundary_values().keys())
                                               
    # Determine the x,y Index in a p-order vector function.
    # Assumption: the function space and vector function space are ordered
    # similarly. Only that the vector function is twice as big because of the 
    # the extra data.
    # For n-D data: :math:`nIndex = iIndices*n + (n-1)`
    xIndices = iIndices*2      # x index is the first data in the 2-D vector function
    yIndices = iIndices*2 + 1  # y-index is the second data in the 2-D vector function
    
    return xIndices, yIndices
    
    
    
def vectorDOF_coordinates(mesh,V,xIndices):
    r"""
    Function to locate the coordinates of the vector functionspace 
    :math:`\mathbf{V} of the mesh **mesh** and indexed according to
    vector DOF boundary index **vectorDOF_boundaryIndex**.
    
    The vector DOF boundary index **vectorDOF_boundaryIndex** can be calculated
    using **vectorDOF_boundaryIndex**.
    
    Usage
    -----
    .. code-block:: python
        
        vectorDOF_coordinates(mesh,V,xIndices)
    
    Parameters
    ----------
    mesh : dolfin.cpp.mesh
           the dolfin mesh of the fluid domain.
           
    V : dolfin.functions.functionspace.VectorFunctionSpace
        the vector function space of the mesh representing a vector-valued 
        finite element function space.
        
    xIndices : numpy.ndarray(float64), shape (N_boundaryDOFs,)
               the :math:`x` indices of the vector boundary DOFs.
               
    Returns
    -------
    xyBoundaryCoordinates : numpy.ndarray(float64), shape (2,N_boundaryDOFs) 
    """
    
    # List all the coordinates of vector function space as ordered by DOFmap
    xyCoordinates = V.dofmap().tabulate_all_coordinates(mesh)

    # Extract vector boundary DOF coordinates    
    xyBoundaryCoordinates = xyCoordinates.reshape(-1,2).T[:,xIndices]
    
    # return vector dof boundary coordinates
    return xyBoundaryCoordinates

                                             

def locate_boundaryDOFs_slow(mesh,boundaryDomains,boundaryID):
    r"""
    
    .... work in progress ....
    
    Compute the mapping from vertices to degrees of freedom only at the 
    boundary. This saves time, since we do not have to loop over all the 
    cells of the domain. If the for loop was replaced by a simple c++ 
    function call, it would be much faster.
    
    
    Usage
    -----
    .. code-block:: python
    
        locate_boundaryDOFs()
    
    
    Parameters
    ----------
    mesh : dolfin.cpp.mesh
           the mesh where we want to compute the mapping between vertices
           and degrees of freedom.
                    
    boundaryDomains : dolfin.cpp.mesh.MeshFunctionSizet
                      A MeshFunction which assigns integer values to the 
                      edges of the mesh.
                                
    exteriodID : int
                 An integer value which we use to identify the 
                 boundaryedges we want to use. The MeshFunction should have
                 this value assigned to some of the edges. 
                                   
    p : int
        the order of the finite element space
        
    M : int
        number of panels. To satisfy, the :math:`M` should be factor of 
        number of body points in the mesh.
        
        
    Assigns
    -------
    boundary_vectorDOF : ndarray (float64), shape (2, nBoundaryDOFs)
                         the array that stores the mapping from vertex 
                         number to degree of freedom number for both scalar 
                         and vector valued functions, of the boundary DOF
                         coordinates.
    
    _boundary_DOFCoordinates_local : ndarray (float64), shape (2, nBoundaryDOFs)
                                     the :math:`x,y` coordinates 
                                     (respectively) of the navier-stokes 
                                     mesh exterior boundary DOFs in the 
                                     *LOCAL* coordinate system.

    _boundary_sortingIndex : ndarray (int64), shape (nBoundaryDOFs,)
                             the array that stores the mapping for sorting
                             the boundary DOF coordinates in anti-clockwise
                             from :math:`-\pi` to :math:`\pi` rad.

    _boundary_DOFCoordinates_local_sorted : ndarray (float64), shape (2, nBoundaryDOFs)
                                            the sorted :math:`x,y` 
                                            coordinates (respectively) of
                                            the navier-stokes mesh exterior 
                                            boundary DOFs in the *LOCAL* 
                                            coordinate system.
    
    _body_DOFCoordinates_local : ndarray (float64), shape (2, nBodyDOFs)
                                 the :math:`x,y` coordinates (respectively)
                                 of the navier-stokes mesh interior boundary
                                 (body) DOFs, in the *LOCAL* coordinate system.
    
    _body_sortingIndex : ndarray (int64), shape (nBodyDOFs,)
                         the array that stores the mapping for sorting the
                         body DOF coordinates in anti-clockwise from 
                         :math:`-\pi` to :math:`\pi` rad.
    
    _body_DOFCoordinates_local_sorted  : ndarray (float64), shape (2, nBodyDOFs)
                                         the sorted :math:`x,y` coordinates
                                         (respectively) of the navier-stokes 
                                         mesh interior boundary (body) DOFs,
                                         in the *LOCAL* coordinate system.
                                        
    
    Returns
    -------
    (-)
    
    
    :First Added:   2013-07-15
    :Last Modified: 2013-12-18
    :Copyright:     Copyright (C) 2013 Artur Palha, Lento Manickathan **pHyFlow**
    :Licence:       GNU GPL version 3 or any later version
                 
    """
    
    # Order of the finite element space
    p = 2
    
    # get the connectivity of the mesh
    #mesh.init(1) # initialize edges
    mesh.init()
    mesh_topology = mesh.topology() # get the mesh topology
    conn12 = mesh_topology(1,2) # get the cells to which an edge belongs
        
    # get the number of edges
    n_edges = boundaryDomains.size()


    # get the indices of the boundary edges, as given by the boundaryFunction
    boundaryEdges = numpy.arange(n_edges)
    boundaryEdges = boundaryEdges[boundaryDomains.array() == boundaryID]
            
    
    # Define the function spaces associated to the mesh
    # for now only CG of order 1 works for sure
    V = dolfin.FunctionSpace(mesh, "CG", p)
    Vv = dolfin.VectorFunctionSpace(mesh, "CG", p, dim=2)

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
    boundary_vectorDOFIndex = numpy.zeros([2,n_dofs], dtype=numpy.uintp)
        
    # Allocate memory space of the arrays that store the coordinates of the dofs
    boundary_DOFCoordinates = numpy.zeros([2,n_dofs])
    
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
        boundary_DOFCoordinates[:,cell_dofs] = dm.tabulate_coordinates(dolfin.Cell(mesh,cell_ind)).T
        
        # now we do the same but for vector valued quantites
        # since, in this case, our vector have two quantities (the x and y components)
        # we have to iterate over the two sub dof maps
        for i, (dms_i, dmcs_i) in enumerate(zip(dms, dms)):
            dms_i.cell_dofs(cell_ind)
            boundary_vectorDOFIndex[i,cell_dofs] = numpy.array(dms_i.cell_dofs(cell_ind),dtype=numpy.uintp)

    # get the list of dofs at the boundary
    scalar_dof_at_boundary = numpy.array(dolfin.DirichletBC(V, dolfin.Constant(0.0), boundaryDomains, boundaryID).get_boundary_values().keys())

    # Assigns vertex to dof mapping
    boundary_DOFCoordinates = boundary_DOFCoordinates[:,scalar_dof_at_boundary]
    boundary_vectorDOFIndex = boundary_vectorDOFIndex[:,scalar_dof_at_boundary]

    # Return coordinates and the its index in function class
    return boundary_DOFCoordinates, boundary_vectorDOFIndex    