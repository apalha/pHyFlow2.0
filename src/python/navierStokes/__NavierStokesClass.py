"""

NAVIERSTOKES
=============

The main class for the navier-stokes method. This class contains of the 
functions related to the calculation of the panel method.

Class structure
---------------
navierSokes:

#. __init__
#. getCoordintes
#. setVelocity
#. setPressure
#. getBoundaryCoordinates
#. evolve
#. getVorticity
#. getMeshPosition
#. getProbe
#. plotVelocity
#. plotPressure
#. plotVorticity
#. save
#. saveVelocity 
#. savePressure
#. saveVorticity

:First Added:   2013-12-10
:Last Modified: 2013-12-27                                                             
:Copyright:     Copyright (C) 2013 Lento Manickathan **pHyFlow**
:License:       GNU GPL version 3 or any later version
"""

#   Copyright (C) 2013 Lento Manickathan                                                                         
#   
#   This file is part of pHyFlow.                                                                                                      
#   
#   pHyFlow is free software: you can redistribute it and/or modify                                                                    
#   it under the terms of the GNU Lesser General Public License as published by                                                       
#   the Free Software Foundation, either version 3 of the License, or                                                                 
#   (at your option) any later version.                                                                                               
#   
#   pHyFlow is distributed in the hope that it will be useful,                                                                         
#   but WITHOUT ANY WARRANTY; without even the implied warranty of                                                                    
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the                                                                      
#   GNU Lesser General Public License for more details.                                                                               
#      
#   You should have received a copy of the GNU Lesser General Public License                                                          
#   along with pHyFlow. If not, see <http://www.gnu.org/licenses/>.    

# External packages
#import dolfin
import fenicstools as _fenicstools
import numpy as _numpy
import copy as _copy

# Import pHyFlow
#from pHyFlow import options
import base as _base

__all__ = ['NavierStokes','solverType','fluidID','noSlipID','extID']

# Default settings: Constants
solverType  = 'ipcs'
fluidID     = 1
noSlipID    = 2
extID       = 3


class NavierStokes(object):
    r"""
    Class containing all the overhead functions and managing the Navier-Stokes 
    solver using FEniCS/DOLFIN Finite Element Solver.
    
    Usage
    -----
    .. code-block:: python
    
        NSDomain = navierStokes(initFile='filePath/fileName.npz')
        
    or
    
    .. code-block:: python

        NSDomain = navierStokes(mesh=mesh, boundaryDomains=boundaryDomains,
                                cmGlobal=cmGlobal, thetaGlobal=thetaGlobal,
                                uMax=uMax, nu=nu, cfl=cfl, 
                                origin=np.array([x0,y0]), L=np.array([Lx,Ly]),
                                N=np.array([nx,ny]))
        
    Parameters
    ----------
    initFile : str
        the filename of the file containing all the parameters to 
        re-initialize the class.

    or
               
    mesh : str
        the mesh data filename.
   
    boundaryDomains : str
        the boundary mesh domain filename.
                  
    cmGlobal : ndarray *(float64)*, shape *(2,)*
        the :math:`x,y` position of the mesh domain in global coordinates.
        
    thetaGlobal : float64
        the rotational angle :math:`\theta` of the mesh domain in global
        coordinates.
        
    uMax : float64
        the maximum fluid velocity :math:`U_{max}`.
        
    nu : float64
        the fluid kinematic viscosity :math:`\nu`.
        
    cfl : float64
        the :math:`CFL` stability parameter.
        
    origin : ndarray *(float64)*, shape *(2,)*
        the :math:`x,y` coordinate of the origin of the probe mesh origin.
        
    L : ndarray *(float64)*, shape *(2,)*
        the width and the height :math:`L_x,L_y` of the probe mesh.
                
    N : ndarray *(float64)*, shape *(2,)*
        the number of probes in the :math:`x,y` direction.
        
    
    Attribute
    ---------
    -
    
    Methods
    -------
    getCoordintes

    setVelocity

    getBoundaryCoordinates
    
    evolve

    getVorticity

    getMeshPosition

    getProbe
    
    Returns
    -------
    -
   
    """
    """
    Revisions
    ---------
    2013-12-10, Lento Manickathan
        - First added


    """
    
    def __init__(self,initFile=None,**kwargs):
        
        
        # Assign Initializing file
        self.__set('initFile', initFile)
        
        # If init file is not given, explicitly define the parameters
        if self.__initFile is None:
            inputData = _copy.deepcopy(kwargs) # Make deep copy
        # If init file is given
        else:
            inputData = _numpy.load(self._initFile) # Load the file
            
        # Initialize the parameters
        self.__set('cmGlobal', inputData['cmGlobal']) # Global reference coordinates
        self.__set('thetaGlobal',inputData['thetaGlobal']) # Local rotation angle
        self.__set('uMax',inputData['uMax']) # Maximum velocity
        self.__set('nu',inputData['nu']) # Kinematic viscosity
        self.__set('cfl',inputData['cfl']) # CFL number
        self.__set('origin',inputData['origin']) # probe mesh grid origin
        self.__set('L',inputData['L']) # probe mesh grid length
        self.__set('N',inputData['N']) # number of probes in x-y direction

        
        # Choose the solver: (Chorin or other solvers)
        # - Various Navier-stokes solving algorithms can be imported. All the
        # algorithms should be stored in './solvers'. The solver is then
        # imported in as part of this module into -> 'self.solver'.
        if solverType == 'chorin':
            self.__solver = _base.chorin(inputData['mesh'],inputData['boundaryDomains'],
                                       self.__nu,self.__cfl,self.__uMax)
        elif solverType == 'ipcs':
            self.__solver = _base.ipcs(inputData['mesh'],inputData['boundaryDomains'],
                                       self.__nu,self.__cfl,self.__uMax)
        else:
            raise NotImplemented('Solver type not implemented or unknown !!')
            
        # Define the probe grid
        self.__xyProbes = _numpy.meshgrid(_numpy.linspace(0,self.__L[0],self.__N[0]) + self.__origin[0],
                                          _numpy.linspace(0,self.__L[1],self.__N[1]) + self.__origin[1])
                     
        # Initialize the probes                     
        self.__probes = _fenicstools.Probes(_numpy.append(self.__xyProbes[0].reshape(-1,1),
                                                          self.__xyProbes[1].reshape(-1,1),axis=1).flatten(),
                                            self.__solver.X)
                                            
            

    def getCoordinates(self):
        r"""
        Function get all the coordinates of the velocity function space
        :math:`\mathbf{V}`. With the returned coordinates, one cold calculate
        the velocity field in the navier-stokes domain.
                
                
        Usage
        -----
        .. code-block:: python
        
            x,y = getCoordinates()
            
        or
        
        .. code-block:: python
        
            xy = getCoordinates()
        
        
        Returns
        -------
        xyCoordinates : ndarray *(float64)*, shape *(2,nDOFs)*
            the :math:`x,y`-coordinates of the vector function space 
            :math:`\mathbf{V}` dofs.
                      
        :First Added:   2013-12-18
        :Last Modified: 2013-12-19
        :Copyright:     Copyright (C) 2013 Lento Manickathan **pHyFlow**
        :Licence:       GNU GPL version 3 or any later version               
        """
        
        # return dof coordinates of the vector function space
        # x,y 1D array
        return self.__solver.V.dofmap().tabulate_all_coordinates(self.__solver.mesh).reshape(-1,2).T[:,::2]
        


    def setVelocity(self,vx,vy):
        r"""
        Function to replace the current `(initial)` velocity Field.
        
        
        Usage
        -----
        .. code-block:: python
            
            setVelocity(vx,vy)
        
        
        Parameters
        ----------
        vx : ndarray *(float64)*, shape *(nDOFs,)*
            the :math:`x`-component of the velocity at the vector function
            space :math:`\mathbf{V}` DOFs.
                     
        vy : ndarray *(float64)*, shape *(nDOFs,)*
            the :math:`y`-component of the velocity at the vector function
            space :math:`\mathbf{V}` DOFs.
        
        
        Assigns
        -------
        vNew : ndarray *(float64)*, shape *(2,nDOFs)*
            the :math:`x,y` velocity at the vector function space 
            :math:`\mathbf{V}` DOFs. 
        
        
        :First Added:   2013-12-18
        :Last Modified: 2013-12-19
        :Copyright:     Copyright (C) 2013 Lento Manickathan **pHyFlow**
        :Licence:       GNU GPL version 3 or any later version        
        """
        
        # vField to assign
        vField = _numpy.vstack((vx,vy)).T.flatten()
        
        # Apply new data
        self.__solver.u1.vector()[:] = vField
            
            
        
    def getBoundaryCoordinates(self):
        """
        Returns the boundary DOF coordinates :math:`x,y_{boundary}` of the 
        vector function space :math:`\mathbf{V}`.
        
        
        Usage
        -----
        .. code-block:: python
        
            xBoundary, yBoundary = getBoundaryCoordinates()
            
        or
        
        .. code-block :: python
        
            xyBoundary = getBoundaryCoordinates()
            
            
        Returns
        -------
        xyBoundary : ndarray *(float64)*, shape *(2,nDOFs)*
            the :math:`x,y` coordinates of the dirichlet boundary where the
            external velocity is applied.
        
        
        :First Added:   2013-12-18
        :Last Modified: 2013-12-18
        :Copyright:     Copyright (C) 2013 Lento Manickathan **pHyFlow**
        :Licence:       GNU GPL version 3 or any later version    
        """
        
        # Get the boundary coordinates of the mesh
        return self.__solver.boundary_DOFCoordinates
        
        

    def evolve(self,vx,vy,cm,theta,cmDot,thetaDot):
        r"""
        Function to evolve the Navier-Stokes by one step with the :math:`x,y`
        velocity boundary conditions at the navier-stokes dirichlet boundary.
        The function will calculate the new velocity and the pressure fields.
        
        To be implemented:
            The new mesh position is used to update the mesh position, whereas
            the current mesh velocity is used to calculate the modified
            convection term to take in account of the rigid mesh motion.
            
        
        Usage
        -----
        .. code-block:: python
        
            evolve(vx,vy,cm,theta,cmDot,thetaDot)
        
        Parameters
        ----------
        vx : ndarray *(float64)*, shape *(nDOFs,)*
            the :math:`x` component of the dirichlet velocity boundary
            condition at the navierstokes DOF boundaries.

        vy : ndarray *(float64)*, shape *(nDOFs,)*
            the :math:`y` component of the dirichlet velocity boundary
            condition at the navierstokes DOF boundaries.                                  
        
        cm : ndarray *(float64)*, shape *(2,)*
            the new :math:`x,y` mesh position of the navierStokes domain.
        
        theta : float64
            the new mesh rotational angle :math:`\theta`.
        
        cmDot : ndarray *(float64)*, shape *(2,)*
            the current :math:`x,y` mesh velocity.
        
        thetaDot : float64
            the current mesh rotational velocity :math:`\dot{\theta}`. 
                   
        Assigns
        -------
        vField : 
            the new velocity field
                
        pField : 
            the new pressure field
        
        Returns
        -------
        -
        
        
        :First Added:   2013-12-19
        :Last Modified: 2013-12-19
        :Copyright:     Copyright (C) 2013 Lento Manickathan **pHyFlow**
        :Licence:       GNU GPL version 3 or any later version   
        """
        
        # Set boundary conditions
        self.__solver.boundaryConditions(vx,vy)
        
        # Solve the navier-stokes problem
        self.__solver.solve()
        


    def getVorticity(self):
        r"""
        Function to evaluat ethe vorticity at the probe coordinates defined by
        the probe mesh
        
        Usage
        -----
        .. code-block:: python
        
            wGrid = getVorticity()
            
        Parameters
        ----------
        -        
        
        Assigns
        -------
        -
        
        Returns
        -------
        wGrid : ndarray *(float64)*, shape *(NxProbes,NyProbes)*
            the vorticity :math:`\omega` at the  probe grid coordinates in the
            navier-stokes domain.
        
        :First Added:   2013-12-22
        :Last Modified: 2013-12-22
        :Copyright:     Copyright (C) 2013 Lento Manickathan **pHyFlow**
        :Licence:       GNU GPL version 3 or any later version   
        """
        
        # Remove old data
        if self.__probes.number_of_evaluations() > 0:
            self.__probes.clear()
        
        # Probe the data
        self.__probes(self._solver.vorticity())
        
        # Return the data
        return self.__probes.array().reshape(self.__N)
        
        
        
    def getProbeGrid(self,metadata=True):
        r"""
        Return the probe grid parameters
        
        Usage
        -----
        .. code-block:: python
            
            origin, L, N = getProbeGrid(metaData=True)
            
        or
        
        .. code-block:: python
        
            xProbes, yProbes = getProbeGrid(metaData=False)
            
        Parameters
        ----------
        -
        
        Assigns
        -------
        -
        
        Returns
        -------
        metadata=True
            origin : ndarray *(float64)*, shape *(2,)*
                the :math:`x,y` coordinate of the origin of the probe mesh
                origin.
             
            L : ndarray *(float64)*, shape *(2,)*
                the width and the height :math:`L_x,L_y` of the probe mesh.
                
            N : ndarray *(float64)*, shape *(2,)*
                the number of probes in the :math:`x,y` direction.
                
        metadata=False
            xProbes : ndarray *(float64)*, shape *(NxProbes,NyProbes)*
                the :math:`x` coordinate of the probe mesh.
                      
            yProbes : ndarray *(float64)*, shape *(NxProbes,NyProbes)*
                the :math:`y` coordinate of the probe mesh.
        
        :First Added:   2013-12-22
        :Last Modified: 2013-12-22
        :Copyright:     Copyright (C) 2013 Lento Manickathan **pHyFlow**
        :Licence:       GNU GPL version 3 or any later version                      
        """
        # Return only the metadata
        if metadata:
            return self.__origin, self.__L, self.__N
        # Return the x,y coordinates
        else:
            return self.__xyProbes
        
        
        
    def __set(self, varName, var):
        """
        Function to check the input parameters
        """
        
        # Check/Set init file
        if varName == 'initFile':
            if type(var) == str or var == None:
                self.__initFile = var
            else:
                raise TypeError('initFile must be a str. It is %s' % str(type(var)))
        
        elif varName == 'cmGlobal':
            if type(var) != _numpy.ndarray:
                raise TypeError('%s must be a numpy ndarray. It is %s.' % (varName, str(type(var))))
            if var.shape != (2,):
                raise ValueError('%s must have shape (%g,0). It has shape %s.' % (varName, 2, str(var.shape)))
            self.__cmGlobal = var
        
        elif varName == 'thetaGlobal':
            if type(var) != float and type(var) != _numpy.float64:
                raise TypeError('%s must be a float. It is %s.' % (varName, str(type(var))))
            self.__thetaGlobal = var
            
        elif varName == 'uMax':
            if type(var) != float and type(var) != _numpy.float64:
                raise TypeError('%s must be a float. It is %s.' % (varName, str(type(var))))
            self.__uMax = var
            
        elif varName == 'nu':
            if type(var) != float and type(var) != _numpy.float64:
                raise TypeError('%s must be a float. It is %s.' % (varName, str(type(var))))
            self.__nu = var
            
        elif varName == 'cfl':
            if type(var) != float and type(var) != _numpy.float64:
                raise TypeError('%s must be a float. It is %s.' % (varName, str(type(var))))
            self.__cfl = var
            
        elif varName == 'origin':
            if type(var) != _numpy.ndarray:
                raise TypeError('%s must be a numpy ndarray. It is %s.' % (varName, str(type(var))))
            if var.shape != (2,):
                raise ValueError('%s must have shape (%g,0). It has shape %s.' % (varName, 2, str(var.shape)))
            self.__origin = var
            
        elif varName == 'L':
            if type(var) != _numpy.ndarray:
                raise TypeError('%s must be a numpy ndarray. It is %s.' % (varName, str(type(var))))
            if var.shape != (2,):
                raise ValueError('%s must have shape (%g,0). It has shape %s.' % (varName, 2, str(var.shape)))
            self.__L = var
            
        elif varName == 'N':
            if type(var) != _numpy.ndarray:
                raise TypeError('%s must be a numpy ndarray. It is %s.' % (varName, str(type(var))))
            if var.shape != (2,):
                raise ValueError('%s must have shape (%g,0). It has shape %s.' % (varName, 2, str(var.shape)))
            self.__N = var

            

    @property
    def cmGlobal(self):
        r"""
        The :math:`x,y` position of the geometry reference point in the global
        coordinate system.
        """                
        return _numpy.copy(self.__cmGlobal)

    @property
    def thetaGlobal(self):
        r"""
        The rotational angle (:math:`rad`) of the geometry in global coordinate
        system.
        """
        return _numpy.copy(self.__thetaGlobal)
        
    @property
    def uMax(self):
        r"""
        The maximum fluid velocity :math:`U_{max}`. It is used to determine the
        maximum time step size :math:`\Delta t_{max}`.
        """
        return _numpy.copy(self.__uMax)
    
    @property
    def nu(self):
        r"""
        The fluid kinematic viscosity :math:`\nu`. For the incompressible ns
        problem, denisty is assumed to be :math:`\rho=1`.
        """
        return _numpy.copy(self.__nu)
        
    @property
    def cfl(self):
        """
        The CFL stability parameter. If explicit time marching scheme, the CFL
        should satisfy :math`CFL < 0`:
        """
        return _numpy.copy(self.__cfl)        
        
    @property
    def origin(self):
        r"""
        The :math:`x,y` coordinate of the origin of the the probe grid.
        """
        return _numpy.copy(self.__origin)

    @property
    def L(self):
        r"""
        The :math:`x,y` size of the probe grid.
        """
        return _numpy.copy(self.__L)

    @property
    def N(self):
        r"""
        The :math:`x,y` number of probe grid nodes.
        """
        return _numpy.copy(self.__N)
        
    @property
    def dtMax(self):
        r"""
        The maximum allowable time step size :math:`\Delta t_{max}` 
        """
        return _numpy.copy(self.__solver.dtMax)
        
        
        