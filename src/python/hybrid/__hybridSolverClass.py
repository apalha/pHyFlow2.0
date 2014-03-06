#-*- coding: utf-8 -*-
__doc__ = """

HYBRID
======

The main class for coupling navier-stokes and vortex-panel class

* Note: No Turbulence scheme implemented (only laminar flow)

Description
-----------
This class wraps the NavierStokes class and VortexPanel class.

Methodology
-----------


Implemented NS Algorithms
-------------------------


:First Added:   2014-02-25
:Last Modified: 2014-02-25                         
:Copyright:     Copyright (C) 2014 Lento Manickathan **pHyFlow**
:License:       GNU GPL version 3 or any later version
"""

#   Copyright (C) 2014 Lento Manickathan                                                                         
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

__all__ = ['HybridSolver']

# External packages
import numpy as _numpy
#import time as _time

# Import point in polygon search function
from matplotlib.nxutils import points_inside_poly as _points_inside_poly
from scipy.interpolate import griddata as _griddata

# Import pHyFlow packages
from pHyFlow.aux.customDecorators import simpleGetProperty
from pHyFlow.blobs import blobOptions
from pHyFlow.lagrangian import LagrangianSolver as _LagrangianSolverClass
from pHyFlow.eulerian import EulerianSolver as _EulerianSolverClass

# Hybrid packages
from pHyFlow.hybrid import hybridOptions as hybridOptions


class HybridSolver(object):
    r"""
    
    Usage
    -----
    .. code-block:: python
    
        Hybrid(lagrangian,multiEulerian,interpolationRegion)
               
        
    Parameters
    ----------
#    blobs : pHyFlow.vortex.VortexBlobs
#            the vortex blob class containing all the parameters that define
#            the vortex-blob problem
#    
#    panels : pHyFlow.panels.Panels
#             the panel class containing all the parameters that define the
#             panel problem.
#    
#    couplingParams : dict, optional
#                     Dictionary containing parameters for coupling the
#                     panels and the blobs together.
#                     
#                     'panelStrengthUpdate' : str
#                                             String parameters specifying if
#                                             panel strength should be updated
#                                             during the mult-step time integration
#                                             scheme (such as RK4).
#                                             
#                                             'constant' : solve for panel strength
#                                                          only in the beginning.
#                                              'varying' : constantly update
#                                                          the panel strength depending
#                                                          on the sub-step position
#                                                          of the particles.
                          
    Attribute
    ---------
#    blobs : pHyFlow.vortex.VortexBlobs
#            the vortex-blobs class is assigned
#            
#    deltaT
#    
#    panels : pHyFlow.panels.Panels
#             the panel blass    
#             
#    t
#    tStep
#    vInf    
#
#    __couplingParams : dict
#                       Dictionary containing parameters for coupling the
#                       panels and the blobs together.  
    
    Methods
    -------
#    evolve
#    evaluateVelocity
#    __coupled_convection
#    __advanceTime
#    __set
#    __mute_blobs_evolve
#    
#    :First Added:   2014-02-26
#    :Last Modified: 2014-02-26                         
#    :Copyright:     Copyright (C) 2014 Lento Manickathan **pHyFlow**
#    :License:       GNU GPL version 3 or any later version   
   
    """
    """
    Revisions
    ---------
    2014-02-26 : Lento Manickathan
                 - First added.

    """
    
    def __init__(self,lagrangian, multiEulerian, interpolationRegions,
                 couplingParams={'adjustLagrangian':hybridOptions.ADJUST_LAGRANGIAN['default'],
                                 'adjustLagrangianAt':hybridOptions.ADJUST_LAGRANGIAN_AT['default'],
                                 'eulerianInitialConditions':hybridOptions.EULERIAN_INITIAL_CONDITIONS['default']},
                 interpolationParams={'algorithm':hybridOptions.INTERPOLATION_ALGORITHM['default'],
                                      'method':hybridOptions.INTERPOLATION_METHOD['default']}):
        
        
        #---------------------------------------------------------------------
        # Check/Set input parameters
            
        # Initialize the parameters
        self.__set_variables('lagrangian',lagrangian) # vortex-blobs class
        self.__set_variables('multiEulerian', multiEulerian) # panel class
        self.__set_variables('interpolationRegions', interpolationRegions)
        self.__set_variables('couplingParams',couplingParams)
        self.__set_variables('interpolationParams',interpolationParams)
        
        #---------------------------------------------------------------------
                
        #---------------------------------------------------------------------
        # Define the initial conditions
        
        # Should we use the solution of the lagrangian domain as the 
        # initial conditions for the eulerian domain
        if self.__couplingParams['eulerianInitialConditions'] == 'lagrangian_field':
            self.__set_eulerianIC_fromLagrangianField()
            
        #---------------------------------------------------------------------            

        #---------------------------------------------------------------------
        # Define the time integration parameters
        
        # Function to modify the eulerian time-step size such that
        # the lagrangian time-step is a multiple of eulerian time step
        self.__set_eulerian_deltaT()
           
        # Get the boundary velocity of the navier-stokes at the initial 
        # time-step.
        #self.__vxyEulerianBoundary = self.navierStokes.getBoundaryVelocity()
        self.__get_eulerian_initialBoundaryVelocity()
        
        #---------------------------------------------------------------------

        #---------------------------------------------------------------------        
        # Verify the coupling parameters

        for i,subDomainID in enumerate(self.multiEulerian.subDomainKeys):
                
            # Make references to probe grid parameters
            probeGridParams = self.multiEulerian[subDomainID].probeGridParams

            # Determine the grid spacing of probe grid
            hGrid = probeGridParams['L']/probeGridParams['N'] # (hx',hy') in local coordinate system

            # Ensure that grid spacing and blob spacing are same
            if _numpy.abs(hGrid[1]-hGrid[0]) > _numpy.spacing(1):
                print "NS sub-domain : %s, hxGrid : %g" % (subDomainID,hGrid[0])
                print "NS sub-domain : %s, hyGrid : %g" % (subDomainID,hGrid[1])
                raise ValueError('Ensure structured grid is uniform and a square !')
                    
            if _numpy.abs(hGrid[0]-self.lagrangian.blobs.h) > _numpy.spacing(1):
                print "NS sub-domain : %s, hGrid : %g, hBlob : %g" % (subDomainID,hGrid[0], self.lagrangian.blobs.h)
                raise ValueError('Ensure grid spacing and blob spacing are the same !')        
        
        #---------------------------------------------------------------------             
            
                
    def evolve(self,cmGlobalNew,thetaLocalNew,cmDotGlobalNew,thetaDotLocalNew):
        r"""
        Function to evolve the coupled hybrid (blobs+panel and eulerian grid)
        problem. 
        
        During the evolution of the vortex blobs, the no-through boundary condition is taken into 
        account by solving the no-slip panel problem.
        
        * Note: moving panels not implemented yet 
        
        Methodology
        -----------

        Algorithm for the hybrid method
          - to couple navier-stokes with vortex-panels
          - see Daeninck _[1] and Stock _[2]
        
        1. 'Correct' the Lagrangian domain.
             'Correct' the strengths of the vortex-particles that are 
             inside the interpolation region of the navier-stokes  domain.
             
             1.1 Interpolate the solution of the navier-stokes domain on
                  to a structured grid (defined by probe mesh). The grid
                  spacing should of nominal particle spacing 
                  :math:`\Delta x_{grid} = \delta_{blob}`. 
         
             1.2 Determine the blobs that are present in the interpolation
                   region. The interpolation region is defined by two
                   polygons: surface polygon and boundary polygon. The 
                   surface polygon is :math:`d_{surf} \Delta x_{grid}` from
                   the surface. The boundary polygon is 
                   :math:`d_{bdry} \Delta x_{grid}` from the exterior eulerian
                   boundary. From Stock: :math:`d_{surf}=3, d_{bdry}=2` for
                   high-Re.
                   
             1.3 'Correct'/'Adjust'/'Reset' the strength of the blobs inside
                    the region.
                   
             1.4 Determine particles that are too-close the body. If
                   the particles are within surface polygon, replace
                   with zero-circulation strength. This comes from the
                   assumption that all the vorticity near the body is
                   captured by the BEM (panels). This ensures that
                   particles does not enter the body during the diffusion
                   process.
            
        2. Advance Lagrangian field (t -> t_{n+1})
             Time step the coupled vortex-panel using a multi-step
             time integration method.
             * Note: initially particles carry zero circulation but
             far-field velocity field is approximated by panels.
             
             2.1 Determine the panel strengths to cancel the slip velocity
                   at the solid boundary (panel collocation points.)
             
             2.2 Determine the aggregate velocity field of vortex blobs,
                   the panels and the free-stream to convect the particles.
                    
             2.3 Convect the blobs. 
                   If panel-strength is to be updated every sub-step, 
                   repeat steps (2.1 and 2.2) in between the RK steps.
                  
             2.4 Diffuse the vortex blobs (when necessary).
             
             2.5 Redistribute and perform population control (when necessary)
              
             2.6 Advance the internal time of vortex-panel class.
             
        3. Determine the Eulerian boundary condition at (t_{n+1}).
              Find the induced velocity acting on the navier-stokes 
              boundary dof coordinates.
              
              3.1 Evaluate the aggregate induced velocity of panels,
                  blobs and free-stream at the navier-stokes dof
                  coordinates.
                  
        4. Advance the eulerian domain (t -> t_{n+1})
             With the boundary conditions at t_{n+1}, solve for the
             velocity field and the pressure field at navier-stokes
             domain.
        
             4.1 Evolve the navier-stokes using the boundary
                  conditions. (update internal time)
                  
        Now, both Eulerian and hybrid domain at :math:`t_{n+1}`
        representing the correct solutions of each sub-domain.
    
        References
        ----------
        .. [1] Daeninck, G. (2006). Development in Hybrid Approaches: 
               Vortex method with known separation location Vortex method 
               with near-wall Eulerian solver RANS-LES coupling. Universit´
               catholique de Louvain.
        .. [2] Stock, M., Gharakhani, A., & Stone, C. (2010). Modeling Rotor
               Wakes with a Hybrid OVERFLOW-Vortex Method on a GPU Cluster. 
               AIAA Applied Aerodynamics …, 1–12. Retrieved from http://arc.aiaa.org/doi/pdf/10.2514/6.2010-4553        
        
           
        Usage
        -----
        .. code-block:: python
        #
        #            evolve()

        Parameters
        ----------
        #        None.
        
        Returned
        --------
        None returned.

        Attributes
        ----------
        #        blobs : pHyFlow.vortex.VortexBlobs
        #                the position and the strength of the blobs are updated
        #                
        #                x : numpy.ndarray(float64), shape (blobs.numBlobs, )
        #                    the :math:`x`-position of the particle.
        #                y : numpy.ndarray(float64), shape (blobs.numBlobs, )
        #                    the :math:`y`-postion of the particles.
        #                g : numpy.ndarray(float64), shape (blobs.numBlobs, )
        #                    the circulation :math:`\Gamma` of the particles.
        #                    
        #        panels : pHyFlow.panels.Panels
        #                 the strength of the panels are updated.
        #                 
        #                 sPanel : numpy.ndarray(float64), shape (panels.nPanelTotal, )
        #                          the strength :math:`\gamma` of the panels.
        #                  
        #                 * Note: position update not implemented yet.
        #         
        #        t : float
        #            the current time
        #            
        #        tStep : int
        #                the current time step
         
        :First Added:   2014-02-26
        :Last Modified: 2014-02-26
        :Copyright:     Copyright (C) 2014 Lento Manickathan, **pHyFlow**
        :Licence:       GNU GPL version 3 or any later version

        """
        
        #----------------------------------------------------------------------
        # Step 1 or 4: 'Correct' the Lagrangian domain.
             
        # (a) Retrieve interpolated vorticity from navier-stokes domain
        # (b) Determine the blobs that are present in the interpolation region
        # (c) 'Correct'/'Adjust'/'Reset' the strength of the blobs inside the region.
        # (d) Determine particles that are too-close the body and 
        #     replace with zero-circulation strength
     
        
        # If lagrangian field is to update at the start 
        if self.__couplingParams['adjustLagrangianAt'] == 'start':
            if self.__couplingParams['adjustLagrangian'] == True:
                self.__correct_blobStrength()   
        #----------------------------------------------------------------------                        

        #----------------------------------------------------------------------            
        # Step 2 : Advance Lagrangian field (t -> t_{n+1})
             
        # (a) Evolve the blobs and panels 
        self.lagrangian.evolve() # TODO: panel position update not implemented !
        #----------------------------------------------------------------------                     
             
        #----------------------------------------------------------------------                         
        # Step 3 : Advance the eulerian domain (t-> t_{n+1})
        
        # The time step size of the eulerian domain is smaller than the
        # time-step size of the lagrangian.        
        # Therefore we advance the navier-stokes eulerian field using
        # k time steps (Forward Euler). The boundary conditions at each 
        # sub-step of the FE is obtained by linear interpolation from the boundary
        # at t_n and t_{n+1}

        # Evolve the navier-stokes        
        self.__multiStep_eulerian(cmGlobalNew,thetaLocalNew,cmDotGlobalNew,thetaDotLocalNew)
        #self.__singleStep_eulerian(cmGlobalNew,thetaLocalNew,cmDotGlobalNew,thetaDotLocalNew)
        

        #----------------------------------------------------------------------        
        # Step 1 or 4 : 'Correct' the Lagrangian domain.

        # If lagrangian field is to update at the end    
        if self.__couplingParams['adjustLagrangianAt'] == 'end':
            if self.__couplingParams['adjustLagrangian'] == True:
                self.__correct_blobStrength()

        #----------------------------------------------------------------------            

        #----------------------------------------------------------------------
        # Step 5 : Additional

        # (a) Check if everything worked, assert t_{ns} == t_{vortex} 
        #        self.vortex.blobs._advanceTime()

        if _numpy.abs(self.lagrangian.t - self.multiEulerian.t) > _numpy.spacing(100):
            raise ValueError('the simulation is not to sync !'\
                             'eulerian t : %g\n '\
                             'lagrangian t : %g' % (self.multiEulerian.t,self.lagrangian.t))
        #----------------------------------------------------------------------

            
    
    def __correct_blobStrength(self):
        r"""

        

        :First Added:   2014-02-26
        :Last Modified: 2014-02-27
        :Copyright:     Copyright (C) 2014 Lento Manickathan, **pHyFlow**
        :Licence:       GNU GPL version 3 or any later version

        """
        # Correct blob strength
        # (b) Determine the blobs that are present in the interpolation region
        # (c) 'Correct'/'Adjust'/'Reset' the strength of the blobs inside the region.
        # (d) Determine particles that are too-close the body and 
        #     replace with zero-circulation strength

        # Algorithm : 
        #   Assumption : blobs are remeshing onto the remeshing grid
        #                before the interpolation
        # 
        #   Step 1: Remove the current set of blobs that are inside the 
        #           interpolation region
        #   Step 2: Generate temporary set of blobs inside interpolation region,
        #           which are in sync with the remeshing grid.
        #   Step 3: Interpolate the vorticity from the navier-stokes interpolation
        #           grid on to the blobs that are inside the interpolation region
        #   Step 4: Added the new set of temporary blobs into the list of 
        #           blobs.

        #----------------------------------------------------------------------
        # Step 1: Removed the current blobs inside the interpolation regions

        # Check the circulation before we remove
        self.__totalCirculationBeforeRemove = self.lagrangian.blobs.g.sum()

        # Remove blobs
        self.__remove_old_blobs() #TODO: conservation of circulation
        
        #----------------------------------------------------------------------
        
        #----------------------------------------------------------------------
        # Step 2: Generate temporary set of blobs inside interpolation
        #         region.        

        # Returns the list of new particles in each sub-domain
        xBlobNewList, yBlobNewList = self.__generate_blob_coordinates()
        #----------------------------------------------------------------------
        
        #----------------------------------------------------------------------
        # Step 3: Interpolate the vorticity from navier-stokes interpolation
        #         grid on to the blobs that are inside the interpolation 
        #         region.

        # Return the list of corrected strengths
        if self.__interpolationParams['algorithm'] == 'scipy_griddata':
            gBlobNewList = self.__scipyInterpolate_strength_from_ns(xBlobNewList, yBlobNewList) #TODO: conservation of circulation
            print 'scipy_griddata !!'
        elif self.__interpolationParams['algorithm'] == 'structured_probes':
            gBlobNewList = self.__interpolate_strength_from_ns(xBlobNewList,yBlobNewList) #TODO: conservation of circulation
            print 'Structured probe interpolation not working!\nNeed to fix the error in structured probe interpolation !!'
            #raise NotImplementedError('Structured probe interpolation not working!\nNeed to fix the error in structured probe interpolation !!')
            gBlobNewList2 = self.__scipyInterpolate_strength_from_nsprobes(xBlobNewList,yBlobNewList) #TODO: conservation of circulation
            
            error = (_numpy.array(gBlobNewList) - _numpy.array(gBlobNewList2))/_numpy.array(gBlobNewList2)
            print error
            print _numpy.abs(error).max()

        elif self.__interpolationParams['algorithm'] == 'scipy_probes':
            print 'scipy interpolation from probes !!'
            gBlobNewList = self.__scipyInterpolate_strength_from_nsprobes(xBlobNewList,yBlobNewList) #TODO: conservation of circulation
            

        #----------------------------------------------------------------------        
        # Determine the change in circulation
        
        # Determine the change in circulation
        self.__epsilon = _numpy.array(self.__totalCirculationRemoved / self.__totalCirculationAdded)

        # Conserve circulation
        #gBlobNewList = _numpy.array(gBlobNewList*self.__epsilon)
        #----------------------------------------------------------------------        

        #----------------------------------------------------------------------
        # Step 4 : Added the new set of temporary blobs into the list of 
        #          blobs.

        # Add the new blobs to lagrangian field        
        self.lagrangian.blobs.addBlobs(_numpy.concatenate(xBlobNewList),
                                   _numpy.concatenate(yBlobNewList),
                                   _numpy.concatenate(gBlobNewList))
        
        # Check the circulation after we add
        self.__totalCirculationAfterAdd = self.lagrangian.blobs.g.sum()
        
        self.__lossInCirculation = self.__totalCirculationAfterAdd - self.__totalCirculationBeforeRemove
        

    def __remove_old_blobs(self):
        r"""
        Function to remove the old blobs inside the interpolation domain.
        
        Description
        -----------
        Step 1 of correcting blob strengths
        
        Methodology
        -----------
        1. Determine the current blobs that are present inside the outer
           boundary of the interpolation region.
           
           1.1 Find blobs that are inside the bounding box of the interpolation
               region.
               
           1.2 Find the blobs (which are inside the interpolation region),
               that are also inside the outer boundary (boundary polygon) of
               the interpolation region.
        
        2. Remove the blobs that are inside the boundary (outer) polygon.
    
        Parameters
        ----------
        
        Returns
        -------
        
        Attributes
        ----------
        
        :First Added:   2014-02-28
        :Last Modified: 2014-02-28
        :Copyright:     Copyright (C) 2014 Lento Manickathan, **pHyFlow**
        :Licence:       GNU GPL version 3 or any later version
        """
        
        # Total circulation removed
        self.__totalCirculationRemoved = []
        
        # Iterate through all the sub domains (each mesh) of the navier-stokes
        for subDomainID in self.multiEulerian.subDomainKeys:
            
            # Make reference to the boundary polygon
            #   Format: closed loop polygon
            #   Shape : (:,2)  [x in column 1, y in column2]
            #   Reasoning : necessary of _points_in_poly function
            xyPolygon = self.__interpolationRegions[subDomainID]['boundaryPolygon']
            
            # Determine the bounding box of the interpolation region
            xyMin_boundingBox = xyPolygon.min(axis=0)
            xyMax_boundingBox = xyPolygon.max(axis=0)
            
            # Make reference to old uncorrected blobs
            xBlob = self.lagrangian.blobs.x
            yBlob = self.lagrangian.blobs.y
            
            # Find the blobs that are inside the bounding box of the 
            # interpolation region
            iBoundingBox = _numpy.where((xBlob > xyMin_boundingBox[0]) & 
                                        (xBlob < xyMax_boundingBox[0]) & 
                                        (yBlob > xyMin_boundingBox[1]) & 
                                        (yBlob < xyMax_boundingBox[1]))[0]
        
            # Determine blobs (which are inside bounding box)
            # that are also inside the interpolation region.                                        
            iInside = _points_inside_poly(_numpy.array([xBlob[iBoundingBox], yBlob[iBoundingBox]]).T,
                                          xyPolygon)

            #TODO: for conservation of circulation
            # Determine the circulation that is removed
            self.__totalCirculationRemoved.append(_numpy.sum(self.lagrangian.blobs.g[iBoundingBox[iInside]]))
        
            if iBoundingBox[iInside].shape[0] !=0:
                # Remove old blobs inside polygon
                self.lagrangian.blobs.removeBlobs(iBoundingBox[iInside])

    
    def __generate_blob_coordinates(self):
        r"""
        Function to generate blobs inside the interpolation region
        
        Description
        -----------
        Step 2 of correcting blob strengths
        
             Generate temporary set of blobs inside interpolation region
            
               (3) +-----------+ (2)
                   |     +     | <------ Bounding Box
                   |   .   .<--|----- interpolation Grid
                   | +       + |         
                   |   .   .   |          y ^
                   |     +     |            |
                   +-----------+           -|---> x
               (0)               (1)
                
        
        
        Methodology
        -----------
        1. Generate blobs inside the bounding box of the interpolation region.
            - The blobs are generated on top of the remeshing grid
            - The blos are generated using mesh grid.
        
        2. Determine which of those blobs are inside the interpolation region.
            - Use the point_in_poly search for determining the which blobs
              are inside both the surface polygon (if it exists) and the 
              boundary polygon
           
        3. Return only the blobs that are inside.
        
        Parameters
        ----------
        
        Returns
        -------
        
        Attributes
        ----------
        
        
        :First Added:   2014-02-28
        :Last Modified: 2014-02-28
        :Copyright:     Copyright (C) 2014 Lento Manickathan, **pHyFlow**
        :Licence:       GNU GPL version 3 or any later version        
        """
        #----------------------------------------------------------------------
        # Determine the parameters of the remeshing grid.
        
        # Reference to blob spacing
        hBlob = self.lagrangian.blobs.h
        
        # Determine the x,y bounds of blob field (to define the remeshing grid)
        xBounds = _numpy.array([-hBlob*1.5, hBlob*1.5])
        yBounds = _numpy.array([-hBlob*1.5, hBlob*1.5])        
        
        #----------------------------------------------------------------------
        
        # List of new particles in each sub-domain of the navier-stokes
        xBlobNewList = []        
        yBlobNewList = []
        
        # Iterate through all the sub domains (each mesh) of the navier-stokes
        for subDomainID in self.multiEulerian.subDomainKeys:        
        
            #----------------------------------------------------------------------
            # Generate mesh-grid of blobs defined by the bounding box        
        
            # Determine the bounding box of the interpolation region
            xyMin_boundingBox = self.__interpolationRegions[subDomainID]['boundaryPolygon'].min(axis=0)
            xyMax_boundingBox = self.__interpolationRegions[subDomainID]['boundaryPolygon'].max(axis=0)        
        
            #----------------------------------------------------------------------
            # Generate temporary set of blobs (in sync with the remeshing grid)
            # 1-D evenly spaced blobs.
            xBlobNew = _numpy.arange(_numpy.floor((xyMin_boundingBox[0] - xBounds[0])/hBlob)-1,
                                     _numpy.ceil( (xyMax_boundingBox[0] - xBounds[0])/hBlob)+1)*hBlob + xBounds[0] + 0.5*hBlob
                                 
            yBlobNew = _numpy.arange(_numpy.floor((xyMin_boundingBox[1] - yBounds[0])/hBlob)-1,
                                     _numpy.ceil( (xyMax_boundingBox[1] - yBounds[0])/hBlob)+1)*hBlob + yBounds[0] + 0.5*hBlob                                 
        
            # Generate 2D grid of blobs
            xBlobNew, yBlobNew = _numpy.meshgrid(xBlobNew, yBlobNew)
        
            # Flatten the blobs
            xBlobNew = xBlobNew.flatten()
            yBlobNew = yBlobNew.flatten()
            #----------------------------------------------------------------------
            
            #----------------------------------------------------------------------
            
            # Concatenate the polygons
            # If there is no surface polygon
            if self.__interpolationRegions[subDomainID]['surfacePolygon'] is None:
                xyPolygon = self.__interpolationRegions[subDomainID]['boundaryPolygon']
            # If there is two polygon lines defining the interpolation region
            else:
                # Concatenate surface and boundary polygon.
                xyPolygon = _numpy.vstack((self.__interpolationRegions[subDomainID]['surfacePolygon'],
                                           self.__interpolationRegions[subDomainID]['boundaryPolygon']))
            
            # Determine which one of them are inside the interpolation grid
            iInside = _points_inside_poly(_numpy.array([xBlobNew, yBlobNew]).T,
                                          xyPolygon)
                                           
            # Append the new particles of the current subdomains to the list 
            # of all the new particles
            xBlobNewList.append(xBlobNew[iInside])
            yBlobNewList.append(yBlobNew[iInside])
        
        # return the list of new particles of all the domains.
        return xBlobNewList, yBlobNewList
        
        
    def __multiStep_eulerian(self,cmGlobalNew,thetaLocalNew,cmDotGlobalNew,thetaDotLocalNew):
        r"""
        Function to advance the navier-stokes eulerian domain.
        
        * Note: navier-stokes motion not implemeneted !
        
        Description
        -----------
        The time step size of the eulerian domain is smaller than the
        time-step size of the lagrangian.        
        Therefore we advance the navier-stokes eulerian field using
        k time steps (Forward Euler). The boundary conditions at each 
        sub-step of the FE is obtained by linear interpolation from the boundary
        at t_n and t_{n+1}
        
        Parameters
        ----------
        
        Returns
        -------
        
        Attributes
        ----------
        
        :First Added:   2014-02-28
        :Last Modified: 2014-02-28
        :Copyright:     Copyright (C) 2014 Lento Manickathan, **pHyFlow**
        :Licence:       GNU GPL version 3 or any later version     
        """
    
        #----------------------------------------------------------------------                         
        # Determine the boundary coordinates of the navier-stokes domains.
        # IMPORTANT : at t_{n+1}
        x,y = self.multiEulerian.getBoundaryCoordinates() # we assume that the coordinates are at t_{n+1}
        #----------------------------------------------------------------------                         
        
        #----------------------------------------------------------------------                         
        # Determine the dirichlet velocity boundary condition at t_{n+1}
        vxNew, vyNew = self.lagrangian.evaluateVelocity(x.copy(),y.copy())
        #----------------------------------------------------------------------
        
        #----------------------------------------------------------------------                         
        # (a) Evolve the navier-stokes problem 
        # The navier-stokes problem is evolved using 'k' sub-steps from t_n to t_{n+1}
        # The boundary condition at each sub-step of the evolve obtained by
        # linear interpolation of the boundary at t_n and t_{n+1}
        # Reference : Daeninck 2006, (3.3.3 Coupling)
        
        # Determine the velocity differential at the boundary
        # We can use this to determine the velocity at each sub steps.
        #dVxyEulerianBoundary = (vxyEulerianBoundaryNew - self.__vxyEulerianBoundary)/self.nEulerianSubSteps
        dVx = (vxNew - self.__vxyEulerianBoundary[0])/self.nEulerianSubSteps
        dVy = (vyNew - self.__vxyEulerianBoundary[1])/self.nEulerianSubSteps
        
        # Multi-step navier-stoke from 1 to nSubSteps
        for k in range(1,self.nEulerianSubSteps+1):
            
            print k
            # Determine the boundary condition at the current sub-step
            #vxyEulerianBoundary_subStep = self.__vxyEulerianBoundary + (dVxyEulerianBoundary * k)
            vx_substep = self.__vxyEulerianBoundary[0] + (dVx * k)
            vy_substep = self.__vxyEulerianBoundary[1] + (dVy * k)
            
            #  Evolve the navier-stokes to the next sub-step
            #self.multiEulerian.evolve(vxyEulerianBoundary_subStep[0], vxyEulerianBoundary_subStep[1],
            #                          cmGlobalNew,thetaLocalNew,cmDotGlobalNew,thetaDotLocalNew) # TODO: moving ns domain not implemented !
            self.multiEulerian.evolve(vx_substep, vy_substep,
                                      cmGlobalNew,thetaLocalNew,cmDotGlobalNew,thetaDotLocalNew) # TODO: moving ns domain not implemented !
                            
        if _numpy.abs(vxNew - vx_substep).max() > _numpy.spacing(1):
            print _numpy.abs(vxNew-vx_substep).max()
            raise ValueError()
        if _numpy.abs(vyNew-vy_substep).max() > _numpy.spacing(1):
            print _numpy.abs(vyNew-vy_substep).max()
            raise ValueError()
            
        # Update the boundary parameters
        #self.__vxEulerianBoundary[] = _numpy.copy(vxyEulerianBoundaryNew) # to be safe
        self.__vxyEulerianBoundary[0] = _numpy.copy(vxNew) # to be safe
        self.__vxyEulerianBoundary[1] = _numpy.copy(vyNew) # to be safe
        #----------------------------------------------------------------------
        

    def __singleStep_navierStokes(self,cmGlobalNew,thetaLocalNew,cmDotGlobalNew,thetaDotLocalNew):
        r"""
        Function to advance the navier-stokes eulerian domain.
        
        * Note: navier-stokes motion not implemeneted !
        
        Description
        -----------
        The time step size of the eulerian domain is smaller than the
        time-step size of the lagrangian.        
        Therefore we advance the navier-stokes eulerian field using
        k time steps (Forward Euler). The boundary conditions at each 
        sub-step of the FE is obtained by linear interpolation from the boundary
        at t_n and t_{n+1}
        
        Parameters
        ----------
        
        Returns
        -------
        
        Attributes
        ----------
        
        :First Added:   2014-02-28
        :Last Modified: 2014-02-28
        :Copyright:     Copyright (C) 2014 Lento Manickathan, **pHyFlow**
        :Licence:       GNU GPL version 3 or any later version     
        """
    
        #----------------------------------------------------------------------                         
        # Determine the boundary coordinates of the navier-stokes domains.
        # IMPORTANT : at t_{n+1}
        x,y = self.multiEulerian.getBoundaryCoordinates() # we assume that the coordinates are at t_{n+1}
        #----------------------------------------------------------------------                         
        
        #----------------------------------------------------------------------                         
        # Determine the dirichlet velocity boundary condition at t_{n+1}
        vxyEulerianBoundary = self.lagrangian.evaluateVelocity(x.copy(),y.copy())
        #----------------------------------------------------------------------                        
        
        #----------------------------------------------------------------------                           
        #  Evolve the navier-stokes to the next step
        self.multiEulerian.evolve(vxyEulerianBoundary[0], vxyEulerianBoundary[1],
                                  cmGlobalNew,thetaLocalNew,cmDotGlobalNew,thetaDotLocalNew) # TODO: moving ns domain not implemented !
                            
        #----------------------------------------------------------------------

    def __set_eulerian_deltaT(self):
        """
        Function to modify the eulerian time-step size.
        
        :First Added:   2014-02-28
        :Last Modified: 2014-02-28
        :Copyright:     Copyright (C) 2014 Lento Manickathan, **pHyFlow**
        :Licence:       GNU GPL version 3 or any later version         
        """        
        
        # Determine the eulerian time-step size
        deltaTE = self.multiEulerian.deltaT
        
        # Determine the lagrangina time-step size
        deltaTL = self.lagrangian.deltaT
        
        # Determine the number of sub-step
        self.nEulerianSubSteps = _numpy.int64(_numpy.ceil(deltaTL/deltaTE))

        # Modify the navier-stokes time-step
        self.multiEulerian.deltaT = deltaTL/self.nEulerianSubSteps
        
        
    def __interpolate_strength_from_ns(self,xBlobNewList,yBlobNewList):
        r"""
        
        #TODO: !!!! BUG HERE !!!! USING SCIPY
        Function to interpolate the circulation (vorticity * area) from the
        navier-stokes domain.
        
        Description
        -----------
        Step 3 of correcting blob strengths
        
        The vorticity field of the navier-stokes is already
        interpolated onto a structured grid. The interpolation of the blob 
        strength from the structured grid is done using a bi-linear interpolation.
        
        The positions of the particles w.r.t to the grid notes are first
        determined. Then the interpolation weight for each particle is determined.
        Using the weights of the bi-linear interpolation the vorticity at
        the blobs are determined.
        
        Finally, the strength of the blobs are equal to the vorticity times
        the grid spacing. This is valid because, we have ensured that the
        grid spacing of the structured grid in the nominal spacing of the 
        vortex blobs.
        
        Algorithm for linear interpolation of vorticity/circulation
        
        Bi-Linear interpolation of vorticity
        
                (4) i,j+1        (3) i+1,j+1
                   +------------+
                   |            |
                   |     _ O    |  O = Location of the blob
                   |     .|     |  + = 4 corners of the linear interpolation
                   |   .        |
                   | .          |
                   +------------+
              (1) i,j            (2) i+1,j
        
        Bi-Linear Interpolation:
        
        .. math::
               \\omega_{blob}\\left(x^{\prime},y^{\\prime}\\right)={\\textstyle \\sum_{i=1}^{4}\
                       {\\omega_{grid}}_i\\cdot h_i\\left(x^{\\prime},y^{\\prime}\\right)},
        
                   where:
                       left( x^\prime-\Delta x \right)\cdot\left( y^\prime-\Delta y \right)}{\Delta x \Delta y},
                       h_2=-\frac{x^\prime\left( y^\prime-\Delta y \right)}{\Delta x \Delta y}
                       h_3=\frac{x^\prime \cdot y^\prime}{\Delta x \Delta y}
                       h_4=-\frac{y^\prime\left( x^\prime-\Delta x \right)}{\Delta x \Delta y}      

        
        Methodology
        -----------
        1. Determine the parameters of the structured grid.
        
        2. 
        
        Parameters
        ----------
        
        Returns
        -------
        
        Attributes
        ----------
        
        :First Added:   2014-02-28
        :Last Modified: 2014-02-28
        :Copyright:     Copyright (C) 2014 Lento Manickathan, **pHyFlow**
        :Licence:       GNU GPL version 3 or any later version      
        """
        
         # Total circulation removed
        self.__totalCirculationAdded = []
        
        # A list of new vortex strengths
        gBlobNewList = []
        
        # Iterate through all the sub domains (each mesh) of the navier-stokes
        for i,subDomainID in enumerate(self.multiEulerian.subDomainKeys):
        
            # Make references to navier-stokes structured grid location
            cmGlobal = self.multiEulerian[subDomainID].cmGlobal
            thetaLocal = self.multiEulerian[subDomainID].thetaLocal
 
            # Make references to probe grid parameters
            probeGridParams = self.multiEulerian[subDomainID].probeGridParams
            
            # Determine the origin in global coordinates
            originGrid = probeGridParams['origin'] + cmGlobal # (x,y) origin in global coordinate system

            
            # Determine the grid spacing of probe grid
            hGrid = probeGridParams['L']/probeGridParams['N'] # (hx',hy') in local coordinate system
            
            
            # Distance of the particles from the grid origin
            LxBlob =   (xBlobNewList[i]-originGrid[0])*_numpy.cos(thetaLocal) \
                     + (yBlobNewList[i]-originGrid[1])*_numpy.sin(thetaLocal)  
                             
            LyBlob = - (xBlobNewList[i]-originGrid[0])*_numpy.sin(thetaLocal) \
                     + (yBlobNewList[i]-originGrid[1])*_numpy.cos(thetaLocal)
            
             # Indexing the particles w.r.t to the origin
            xLIndex = _numpy.int64(_numpy.floor(LxBlob/hGrid[0]))
            yLIndex = _numpy.int64(_numpy.floor(LyBlob/hGrid[1]))
                

            # Distance of each blob from it's lower left grid box
            xPrime = LxBlob - xLIndex*hGrid[0]# - 0.5*hGrid[0]
            yPrime = LyBlob - yLIndex*hGrid[1]# - 0.5*hGrid[0]

            # Linear Interpolation weights
            h1 =    (xPrime - hGrid[0]) * (yPrime-hGrid[1]) / (hGrid[0]*hGrid[1])
            h2 =       - xPrime         * (yPrime-hGrid[1]) / (hGrid[0]*hGrid[1])
            h3 =         xPrime         *     yPrime        / (hGrid[0]*hGrid[1])
            h4 =       - yPrime         * (xPrime-hGrid[0]) / (hGrid[0]*hGrid[1])      
                
            # Retrieve the interpolated vorticity from ns grid to interpolation grid
            wGrid = _numpy.copy(self.multiEulerian[subDomainID].getVorticity())
                
            # Interpolating the (vorticity * grid area) = circulation from the grid the blobs       
            gBlobNew = (h1*wGrid[xLIndex,yLIndex]       +\
                        h2*wGrid[xLIndex,yLIndex+1]     +\
                        h3*wGrid[xLIndex+1,yLIndex+1]   +\
                        h4*wGrid[xLIndex+1,yLIndex]) * (self.lagrangian.blobs.h**2) # Circulation = vort * area     

            gBlobNewList.append(gBlobNew)
            
        self.__totalCirculationAdded = _numpy.sum(gBlobNewList, axis=1)
            
        return gBlobNewList


    def __scipyInterpolate_strength_from_nsprobes(self, xBlobNewList, yBlobNewList):
        r"""
        Use the scipy interpolate function to interpolate from probe grid
        
        
        :First Added:   2014-03-05
        :Last Modified: 2014-03-05
        :Copyright:     Copyright (C) 2014 Lento Manickathan, **pHyFlow**
        :Licence:       GNU GPL version 3 or any later version      
        """
        
        # Total circulation removed
        self.__totalCirculationAdded = []
        
        # A list of new vortex strengths
        gBlobNewList = []
        
        # Iterate through all the sub domains (each mesh) of the navier-stokes
        for i,subDomainID in enumerate(self.multiEulerian.subDomainKeys):
        
            # Stack coordinates
            x,y = self.multiEulerian.subDomains[subDomainID].probeGridMesh
            
            # Reshape the data to (M,2)
            xy = _numpy.vstack((x.flatten(),y.flatten())).T
            
            # Retrieve the vorticity
            w = self.multiEulerian.subDomains[subDomainID].getVorticity()
            
            # Reshape the vorticity data
            ww = w.flatten()
            
            # Interpolate vorticity to grid
            wBlobNew = _griddata(xy,ww,(xBlobNewList[i],yBlobNewList[i]), method=self.__interpolationParams['method'])
                 
            gBlobNew = wBlobNew * (self.lagrangian.blobs.h**2)
            
            gBlobNewList.append(gBlobNew)
            
        self.__totalCirculationAdded = _numpy.sum(gBlobNewList, axis=1)
            
        return gBlobNewList  
        
        
    def __scipyInterpolate_strength_from_ns(self,xBlobNewList,yBlobNewList):
        r"""
        Function to interpolate the circulation (vorticity * area) from the
        navier-stokes domain using SCIPY GRIDDATA function.
        
        Description
        -----------
        Step 3 of correcting blob strengths
        
  
        
        Methodology
        -----------
        1. Determine the parameters of the structured grid.
        
        2. 
        
        Parameters
        ----------
        
        Returns
        -------
        
        Attributes
        ----------
        
        :First Added:   2014-03-03
        :Last Modified: 2014-03-03
        :Copyright:     Copyright (C) 2014 Lento Manickathan, **pHyFlow**
        :Licence:       GNU GPL version 3 or any later version      
        """
        
         # Total circulation removed
        self.__totalCirculationAdded = []
        
        # A list of new vortex strengths
        gBlobNewList = []
        
        # Iterate through all the sub domains (each mesh) of the navier-stokes
        for i,subDomainID in enumerate(self.multiEulerian.subDomainKeys):
        
            # Reference the solver
            solver = self.multiEulerian.subDomains[subDomainID]._EulerianSolver__solver
            
            # Extract the coordinates
            x,y = solver.X.dofmap().tabulate_all_coordinates(solver.mesh).reshape(-1,2).T.copy()
            
            # Reshape
            xy = _numpy.hstack((x.reshape(-1,1),y.reshape(-1,1)))
            
            # Extract vorticity
            w = solver.vorticity().vector().array()
            
            # Interpolate vorticity to grid
            wBlobNew = _griddata(xy,w,(xBlobNewList[i],yBlobNewList[i]), method=self.__interpolationParams['method'])
                 
            gBlobNew = wBlobNew * (self.lagrangian.blobs.h**2)
            
            gBlobNewList.append(gBlobNew)
            
        self.__totalCirculationAdded = _numpy.sum(gBlobNewList, axis=1)
            
        return gBlobNewList        
        
    
    def __get_eulerian_initialBoundaryVelocity(self):
        r"""        
        Function to determine the initial eulerian boundary conditions.
        This is vital for the multi-step eulerian scheme.
        """
        
        #----------------------------------------------------------------------
        # Determine the boundary coordinates of the navier-stokes.
        x,y = self.multiEulerian.getBoundaryCoordinates() # we assume that the coordinates are at t_{n+1}
        
        #----------------------------------------------------------------------
        # Determine the induced velocity at the navier-stokes boundary at t_{n}
        vx, vy = self.lagrangian.evaluateVelocity(x.copy(),y.copy())
        
        #----------------------------------------------------------------------
        # Allocate variable space
        self.__vxyEulerianBoundary = _numpy.zeros((2,x.shape[0]))
        
            
        # Store the induced velocity
        self.__vxyEulerianBoundary[0] = _numpy.copy(vx)
        self.__vxyEulerianBoundary[1] = _numpy.copy(vy)

        
    def __set_eulerianIC_fromLagrangianField(self):
        r"""
        Function to replace the eulerian (navier-stokes) initial conditions
        with the lagrangian solutions.
        
        Description
        -----------
        The velocity field from the vortex blobs, panels and free-stream is
        determined at the vector function space DOF coordinates of the 
        navier-stokes grid.
        
        This velocity field is then assigned to the navier-stokes domain.
        
        Methodolgy
        ----------
        1. Determine the vector DOF coordinates of the navier-stokes grids.
        
        2. Determine the induced velocity (blobs + panels + free-stream)
            on the these DOFs.
            
        3. Replace the velocity field of the navier-stokes with the new 
            velocity field.
            
        Parameters
        ----------
        
        Returns
        -------
        
        Assigns
        -------
        
        :First Added:   2014-02-28
        :Last Modified: 2014-02-28
        :Copyright:     Copyright (C) 2014 Lento Manickathan, **pHyFlow**
        :Licence:       GNU GPL version 3 or any later version      
        """
        
        # Determine the vector DOF coordiantes of navier-stokes grid
        x, y = self.multiEulerian.getCoordinates()
        
        # Determine the induced velocity of vortex on these points
        vx, vy  = self.lagrangian.evaluateVelocity(x.copy(),y.copy())
        
        # Replace the origin velocity field.
        self.multiEulerian.setVelocity(vx,vy)
        
        
    def __set_variables(self,varName,var):
        """
        Function to set all the parameters.
        
        Usage
        -----
        .. code-block :: python
        
            __set('varName', var)
            
        Parameters
        ----------
        varName : str
                  the variable name that is to be assigned. Parameters that
                  are set:

        var : <var>
              the variable that is to be set.

        Returns              
        -------
        None returned.
        
        Attributes
        ----------
        #        blobs : pHyFlow.vortex.VortexBlobs
        #                the vortex-blobs class is assigned
        #                
        #        panels : pHyFlow.panels.Panels
        #                 the panel blass
        #                 
        #        couplingParams : dict
        #                         Dictionary containing parameters for coupling the
        #                         panels and the blobs together.
                         
        :First Added:   2014-02-24
        :Last Modified: 2014-02-24
        :Copyright:     Copyright (C) 2014 Lento Manickathan, **pHyFlow**
        :Licence:       GNU GPL version 3 or any later version    
        """
        #----------------------------------------------------------------------
        # navier-stokes
        if varName == 'multiEulerian':
            if str(type(var)) != "<class '__main__.MultiEulerianSolver'>":
                raise ValueError("'navierStokes' should be of type pHyFlow.navierStokes.NavierStokes. It is %s" % type(var))
                
            self.multiEulerian = var
            
        # vortex-panels
        elif varName == 'lagrangian':
            if type(var) != _LagrangianSolverClass:
                raise ValueError("'vortex' should be of type pHyFlow.vortexPanel.VortexPanel."\
                                 "It is %s" % type(var))
            # Everything good 
            self.lagrangian = var 
        #----------------------------------------------------------------------

        #---------------------------------------------------------------------
        # Check and set interpolationRegion parameters
        elif varName == 'interpolationRegions':
                        
            # Check the data type            
            if type(var) != dict:
                raise TypeError('interpolationRegions must be dict of navier-stokes'\
                                ' domains containing dictionary of the surface polygon'\
                                ' and boundary polygon. It is %s.' % type(var))
                                
            # Check if the dictionary contains the interpolation polygon of
            # navier-stokes domains
            for key in self.multiEulerian.subDomainKeys:
                if key not in var.keys():
                    raise ValueError('interpolationRegions[%s] cannot be found.' % key)
            
                # Check the data type            
                if type(var[key]) != dict:
                    raise TypeError('interpolationRegions[%s] must be dict, containing'\
                                    ' surfacePolygon and boundaryPolygon. It is %s.' % type(var))
            
                # Check if surface polygon is described
                if 'surfacePolygon' not in var[key]:
                    raise ValueError("interpolationRegion[%s]['surfacePolygon'] cannot be found" % key)                    
                    
                # Check data of surfacePolygon
                if var[key]['surfacePolygon'] is not None:
                    if type(var[key]['surfacePolygon']) != _numpy.ndarray or var[key]['surfacePolygon'].dtype != float:
                        raise ValueError("interpolationRegions[%s]['surfacePolygon']"\
                                         " should be numpy ndarray of float64. It is %s." % (key,type(var[key]['surfacePolygon'])))
                    
                    # Check data shape
                    if var[key]['surfacePolygon'].shape[1] != 2:
                        raise ValueError("interpolationRegions[%s]['surfacePolygon'] must be of shape (nPoints,2)."\
                                         "It has shape %s." % (key, str(var[key]['surfacePolygon'].shape)))
   
                # Check if surface polygon is described
                if 'boundaryPolygon' not in var[key]:
                    raise ValueError("interpolationRegion[%s]['boundaryPolygon'] cannot be found" % key)                    
                    
                # Check data of surfacePolygon
                if type(var[key]['boundaryPolygon']) != _numpy.ndarray or var[key]['boundaryPolygon'].dtype != float:
                    raise ValueError("interpolationRegions[%s]['boundaryPolygon']"\
                                     " should be numpy ndarray of float64. It is %s." % (key,type(var[key]['boundaryPolygon'])))
                    
                    # Check data shape
                    if var[key]['boundaryPolygon'].shape[1] != 2:
                        raise ValueError("interpolationRegions[%s]['boundaryPolygon'] must be of shape (nPoints,2)."\
                                         "It has shape %s." % (key, str(var[key]['boundaryPolygon'].shape)))                                      
                                         
        
               
            # Everything okay
            self.__interpolationRegions = var
        #----------------------------------------------------------------------

        #----------------------------------------------------------------------            
        # coupling parameters
        elif varName == 'couplingParams':
            if type(var) != dict:
                raise ValueError("'couplingParams' should be of type dict. It is %s" % type(var))

            if var['adjustLagrangian'] not in hybridOptions.ADJUST_LAGRANGIAN['available']:
                raise ValueError("The given couplingParams['adjustLagrangian'] is unknown."\
                                 " It should be in [%s]. It is %s"\
                                 % (str(hybridOptions.ADJUST_LAGRANGIAN['available']),var['adjustLagrangian']))

            if var['adjustLagrangianAt'] not in hybridOptions.ADJUST_LAGRANGIAN_AT['available']:
                raise ValueError("The given couplingParams['adjustLagrangianAt'] is unknown."\
                                 " It should be in [%s]. It is %s"\
                                 % (str(hybridOptions.ADJUST_LAGRANGIAN_AT['available']),var['adjustLagrangianAt']))


            if var['eulerianInitialConditions'] not in hybridOptions.EULERIAN_INITIAL_CONDITIONS['available']:
                raise ValueError("The given couplingParams['eulerianInitialConditions'] is unknown."\
                                 " It should be in [%s]. It is %s"\
                                 % (str(hybridOptions.EULERIAN_INITIAL_CONDITIONS['available']),var['eulerianInitialConditions']))
                                 
            self.__couplingParams = var
                                
        #---------------------------------------------------------------------
        elif varName == 'interpolationParams':
            if type(var) != dict:
                raise ValueError("'interpolationParams' should be of type dict. It is %s" % type(var))

            if var['algorithm'] not in hybridOptions.INTERPOLATION_ALGORITHM['available']:
                raise ValueError("The given interpolationParams['algorithm'] is unknown."\
                                 " It should be in [%s]. It is %s"\
                                 % (str(hybridOptions.INTERPOLATION_ALGORITHM['available']),var['algorithm']))

            if var['method'] not in hybridOptions.INTERPOLATION_METHOD['available']:
                raise ValueError("The given interpolationParams['method'] is unknown."\
                                 " It should be in [%s]. It is %s"\
                                 % (str(hybridOptions.INTERPOLATION_METHOD['available']),var['method']))  
                                 
            self.__interpolationParams = var
            
    #--------------------------------------------------------------------------
    # All attributes
            
    #    @simpleGetProperty
    #    def tStep(self):
    #        r"""
    #        tStep : int
    #                the current time step of the simulation
    #        """
    #        return self.blobs.tStep
    #            
    #    # Current time t
    #    @simpleGetProperty
    #    def t(self):
    #        r"""
    #        t : float
    #            the current time of the simulation
    #        """
    #        return self.blobs.t
    #    
    #    # Time-step size
    #    @simpleGetProperty
    #    def deltaT(self):
    #        r"""
    #        deltaT : float
    #                 the time step size of the simulation :math:`|Delta t`
    #        """
    #        return self.blobs.deltaTc
    #    
    #    # Free-stream velocity
    #    @simpleGetProperty
    #    def vInf(self):
    #        """
    #        vInf : numpy.ndarray(float64), shape (2,)
    #               the :math:`x,y` component of the free-stream velocity
    #        """
    #        return self.blobs.vInf
        
      
       
#    def evaluateVelocity(self,xTarget,yTarget):
#        """
#        Function to evaluate the total induced velocity due to vortex blobs,
#        panels and the free-stream flow.
#        .. math::
#        
#            \mathbf{u} = \mathbf{u}_{\omega} + \mathbf{u}_{\gamma} + \mathbf{u}_\infty
#            
#        Usage
#        -----
#        .. code-block :: python
#        
#            vx, vy = evaluteVelocity(xTarget, yTarget)
#            
#        Parameters
#        ----------
#        xTarget : numpy.ndarray(float64), shape (nTargets,)
#                  the :math:`x`-coordinate of the target location, where the
#                  total velocity is to be evaluated
#                  
#        yTarget : numpy.ndarray(float64), shape (nTargets,)
#                  the :math:`y`-coordinate of the target location, where the
#                  total velocity is to be evaluated
#          
#        Returns
#        -------
#        vx : numpy.ndarray(float64), shape (nTarget,)
#             the :math:`x`-component of the total induced velocity on each
#             **(xTarget, yTarget)** point.
#
#        vy : numpy.ndarray(float64), shape (nTarget,)
#             the :math:`y`-component of the total induced velocity on each
#             **(xTarget, yTarget)** point.
#
#        Attributes
#        ----------
#        None changed.
#        
#        :First Added:   2014-02-24
#        :Last Modified: 2014-02-24                         
#        :Copyright:     Copyright (C) 2014 Lento Manickathan **pHyFlow**
#        :License:       GNU GPL version 3 or any later version   
#        """
#        # Determine the induced velocity due to blobs + panels + free-stream.
#        vx,vy = self.blobs.evaluateVelocity(xTarget,yTarget) \
#                            + self.panels.evaluateVelocity(xTarget,yTarget) \
#                            + self.vInf.reshape(2,-1)
#        
#        # return the induced velocity
#        return vx,vy      