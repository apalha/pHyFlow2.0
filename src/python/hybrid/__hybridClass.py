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

__all__ = ['Hybrid']

# External packages
import numpy as _numpy
from matplotlib.nxutils import points_inside_poly as _points_inside_poly
# Import pHyFlow packages
#from pHyFlow import options as _pHyFlowOptions
#from pHyFlow.aux.customDecorators import simpleGetProperty

# Hybrid packages
#from pHyFlow.hybrid import hybridOptions as _hybridOptions

# Other packages
from pHyFlow.vortexPanel import VortexPanel as _VortexPanel
from pHyFlow.navierStokes import NavierStokes as _NavierStokes


class Hybrid(object):
    r"""
    
    Usage
    -----
    .. code-block:: python
    
        Hybrid(vortexPanel,NavierStokes,Motion,
               )
               
        
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
    
    def __init__(self,vortexPanel, navierStokes, interpolationRegion):
        
        
        #---------------------------------------------------------------------
        # Check/Set input parameters
            
        # Initialize the parameters
        self.__set('vortexPanel',vortexPanel) # vortex-blobs class
        self.__set('navierStokes', navierStokes) # panel class
        self.__set('interpolationRegion', interpolationRegion)            
        
        
                
    def evolve(self,cmGlobalNew,thetaLocalNew,cmDotGlobalNew,thetaDotLocalNew):
        r"""
        Function to evolve the coupled hybrid (vortex-panel and navierStokes)
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
        # Step 1 : 'Correct' the Lagrangian domain.
             
        # (a) Retrieve interpolated vorticity from navier-stokes domain
        # (b) Determine the blobs that are present in the interpolation region
        # (c) 'Correct'/'Adjust'/'Reset' the strength of the blobs inside the region.
        # (d) Determine particles that are too-close the body and 
        #     replace with zero-circulation strength
        self.__correct_blobStrength()    
        #----------------------------------------------------------------------                        

        #----------------------------------------------------------------------            
        # Step 2 : Advance Lagrangian field (t -> t_{n+1})
             
        # (a) Evolve the blobs.
             
        #----------------------------------------------------------------------                     
             
        #----------------------------------------------------------------------                         
        # Step 3 : Determine the Eulerian boundary condition at (t_{n+1}).
        
        # (a) Find the aggregate induced velocity of panels, blobs and 
        #     free-stream at the navier-stokes dof coordinates.
        #x,y = multiNavierstokes.getBoundaryCoordinates()
        #----------------------------------------------------------------------
        # Step 4 : Advance the eulerian domain (t -> t_{n+1})

        # (a) Evolve the navier-stokes problem          
        
        #----------------------------------------------------------------------
        
        #----------------------------------------------------------------------
        # Step 5 : Additional

        # (a) Check if everything worked, assert t
        
        #----------------------------------------------------------------------
        pass
    
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

        #        #----------------------------------------------------------------------
        #        # (a) Retrieve interpolated vorticity from navier-stokes domain                
        #        wGrid = self.navierStokes.getVorticity()
        #
        #        #----------------------------------------------------------------------        

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
        # Define the interpolation region:
        #   The interpolation region polygon and it's bounding box

        # Interpolation region (concatenate surface polygon and boundary polygon)
        #   Format: closed loop polygons, stacked on top of each other
        #   Shape : (:,2)  [x in column 1, y in column2]
        #   Reasoning : necessary of _points_in_poly function
        xyInterpPolygon = _numpy.vstack((self.navierStokes.surfacePolygon,  #TODO: change to global
                                         self.navierStokes.boundaryPolygon)) #TODO: change to global

        # Determine the limit of the interpolation region
        #   Determine the global coordinates of the bounding box of the interp. region
        xyMin_interpRegion = xyInterpPolygon.min(axis=0)
        xyMax_interpRegion = xyInterpPolygon.min(axis=0)

        #----------------------------------------------------------------------

        #----------------------------------------------------------------------
        # Remove old blobs

        # Make reference to old uncorrected blobs
        xBlob, yBlob = self.vortexPanel.blobs.x, self.vortexPanel.blobs.y
        
        # Determine the old blobs inside the boundary box of interpolation region
        iBoundingBox = _numpy.where((xBlob > xyMin_interpRegion[0]) & 
                                    (xBlob < xyMax_interpRegion[0]) & 
                                    (yBlob > xyMin_interpRegion[1]) & 
                                    (yBlob < xyMax_interpRegion[1]))[0]
                                        
        # Determine blobs (which are inside bounding box), that are also inside the interpolation region.                                        
        iInside = _points_inside_poly(_numpy.array([xBlob[iBoundingBox], yBlob[iBoundingBox]]).T,
                                      xyInterpPolygon)

        #TODO: for conservation of circulation
        # Determine the circulation that is removed
        totalCirculationRemoved = _numpy.sum(self.vortexPanel.blobs.g[iBoundingBox[iInside]])
        
        # Remove old blobs inside polygon
        self.vortexPanel.blobs.removeBlobs(self,iBoundingBox[iInside])                                           

        # clear unnecessary parameters
        del iBoundingBox
        del iInside  
        
        #----------------------------------------------------------------------
        
        #----------------------------------------------------------------------
        # Generate temporary set of blobs inside interpolation region
        #
        #   (3) +-----------+ (2)
        #       |     +     | <------ Bounding Box
        #       |   .   .<--|----- interpolation Grid
        #       | +       + |
        #       |   .   .   |
        #       |     +     |  
        #       +-----------+
        #   (0)               (1)
        #        
        
        # Reference to blob spacing
        hBlob = self.vortexPanel.blobs.h
        
        # Determine the x,y bounds of blob field (to define the remeshing grid)
        xBounds = _numpy.array([-hBlob*1.5, hBlob*1.5])
        yBounds = _numpy.array([-hBlob*1.5, hBlob*1.5])        
        
        # Generate temporary set of blobs (in sync with the remeshing grid)
        # 1-D evenly spaced blobs.
        xBlobNew = _numpy.arange(_numpy.floor((xyMin_interpRegion[0] - xBounds[0])/hBlob)-1,
                                 _numpy.ceil( (xyMax_interpRegion[0] - xBounds[0])/hBlob)+1)*hBlob + xBounds[0] + 0.5*hBlob
                                 
        yBlobNew = _numpy.arange(_numpy.floor((xyMin_interpRegion[1] - yBounds[0])/hBlob)-1,
                                 _numpy.ceil( (xyMax_interpRegion[1] - yBounds[0])/hBlob)+1)*hBlob + yBounds[0] + 0.5*hBlob                                 
        
        # Generate 2D grid of blobs
        xBlobNew, yBlobNew = _numpy.meshgrid(xBlobNew, yBlobNew)
        
        # Flatten the blobs
        xBlobNew = xBlobNew.flatten()
        yBlobNew = yBlobNew.flatten()
        
        # Determine which one of them are inside the interpolation grid
        iInside = _points_inside_poly(_numpy.array([xBlobNew, yBlobNew]).T,
                                      xyInterpPolygon)
                                           
        # New particles that are inside polygon
        xBlobNew = xBlobNew[iInside]
        yBlobNew = yBlobNew[iInside]
        
        # Delete unnecessary variables (no need for interpolation region)
        del xyInterpPolygon
        del xyMin_interpRegion
        del xyMax_interpRegion
        del xBounds
        del yBounds
        del iInside

        #----------------------------------------------------------------------
        
        #----------------------------------------------------------------------
        # Interpolate the vorticity from the navier-stokes interpolation grid
        
        # Algorithm for linear interpolation of vorticity
        # Bi-Linear interpolation of vorticity
        #
        #        (4) i,j+1        (3) i+1,j+1
        #           +------------+
        #           |            |
        #           |     _ O    |  O = Location of the blob
        #           |     .|     |  + = 4 corners of the linear interpolation
        #           |   .        |
        #           | .          |
        #           +------------+
        #      (1) i,j            (2) i+1,j
        #
        #   Bi-Linear Interpolation:
        #   .. math::
        #       \\omega_{blob}\\left(x^{\prime},y^{\\prime}\\right)={\\textstyle \\sum_{i=1}^{4}\
        #               {\\omega_{grid}}_i\\cdot h_i\\left(x^{\\prime},y^{\\prime}\\right)},
        #
        #           where:
        #               left( x^\prime-\Delta x \right)\cdot\left( y^\prime-\Delta y \right)}{\Delta x \Delta y},
        #               h_2=-\frac{x^\prime\left( y^\prime-\Delta y \right)}{\Delta x \Delta y}
        #               h_3=\frac{x^\prime \cdot y^\prime}{\Delta x \Delta y}
        #               h_4=-\frac{y^\prime\left( x^\prime-\Delta x \right)}{\Delta x \Delta y}       

        # Make references to navier-stokes grid location
        cmGlobal = self.navierStokes.cmGlobal
        thetaLocal = self.navierStokes.thetaLocal

        # Make references to probe grid parameters
        probeGridParams = self.navierStokes.probeGridParams
        
        # Determine the origin in global coordinates
        originGrid = probeGridParams['origin'] + cmGlobal # (x,y) origin in global coordinate system
        
        # Determine the grid spacing of probe grid
        hGrid = probeGridParams['L']/probeGridParams['N'] # (hx',hy') in local coordinate system
        
        
        # Distance of the particles from the grid origin
        LxBlob =   (xBlobNew-originGrid[0])*_numpy.cos(thetaLocal) \
                 + (yBlobNew-originGrid[1])*_numpy.sin(thetaLocal)
                     
        LyBlob = - (xBlobNew-originGrid[0])*_numpy.sin(thetaLocal) \
                 + (yBlobNew-originGrid[1])*_numpy.cos(thetaLocal)
                 
        # Indexing the particles w.r.t to the origin
        xLIndex = _numpy.int64(_numpy.floor(LxBlob/hGrid[0]))
        yLIndex = _numpy.int64(_numpy.floor(LyBlob/hGrid[1]))                 
        
        # Distance of each blob from it's lower left grid box
        xPrime = LxBlob - xLIndex*hGrid[0]
        yPrime = LyBlob - yLIndex*hGrid[1]        
        
        # Linear Interpolation weights
        h1 =    (xPrime - hGrid[0]) * (yPrime-hGrid[1]) / (hGrid[0]*hGrid[1])
        h2 =       - xPrime         * (yPrime-hGrid[1]) / (hGrid[0]*hGrid[1])
        h3 =         xPrime         *     yPrime        / (hGrid[0]*hGrid[1])
        h4 =       - yPrime         * (xPrime-hGrid[0]) / (hGrid[0]*hGrid[1])        
        
        # Retrieve the interpolated vorticity from ns grid to interpolation grid
        wGrid = self.navierStokes.getVorticity()
        
        # Interpolating the (vorticity * grid area) = circulation from the grid the blobs
        # * Note: we must follow Stock's definition of the grid spacing
        #       grid spacing = nominal particle spacing
        #           :math:`\Delta x_{grid} = \sigma_{blob}` 
        gBlobNew = (h1*wGrid[yLIndex,xLIndex]       +\
                    h2*wGrid[yLIndex,xLIndex+1]     +\
                    h3*wGrid[yLIndex+1,xLIndex+1]   +\
                    h4*wGrid[yLIndex+1,xLIndex]) * (hGrid[0]*hGrid[1]) # Circulation = vort * area        
        
        # Total circulation that is to be added
        totalCirculationAdded = _numpy.sum(gBlobNew)
        
        # Remove unnecessary variables
        del originGrid
        del hGrid
        del LxBlob
        del LyBlob
        del xLIndex
        del yLIndex
        del xPrime
        del yPrime 
        del h1
        del h2
        del h3
        del h4
        
        #----------------------------------------------------------------------
        
        #----------------------------------------------------------------------
        # TODO:  Ensure the conservation of circulation
        
        # Check the change in circulation
        dGamma = totalCirculationRemoved - totalCirculationAdded
        print dGamma
        
        # Ensure circulation does not increase/decrease.
        #gBlobNew = gBlobNew + (dGamma * _numpy.abs(gBlobNew) / _numpy.abs(gBlobNew).sum())
        
        #----------------------------------------------------------------------
        
        #----------------------------------------------------------------------
        # Added the new blobs to list of all blobs
        
        self.vortexPanel.blobs.addBlobs(xBlobNew,yBlobNew,gBlobNew)
        #----------------------------------------------------------------------
        
        
    def __set(self,varName,var):
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
        # navier-stokes
        if varName == 'navierStokes':
            if type(var) != _NavierStokes:
                raise ValueError("'navierStokes' should be of type pHyFlow.navierStokes.NavierStokes. It is %s" % type(var))
                
            self.navierStokes = var
            
        # vortex-panels
        elif varName == 'vortexPanel':
            if type(var) != _VortexPanel:
                raise ValueError("'vortexPanel' should be of type pHyFlow.vortexPanel.VortexPanel. It is %s" % type(var))
            
            self.vortexPanel = var 

        #---------------------------------------------------------------------
        # Check and set interpolationRegion parameters
        elif varName == 'interpolationDomain':
                        
            # Keys in interpolationDomain
            interpolationDomainKeys = ['surfacePolygon','boundaryPolygon']
            
            # check the interpolationDomain data type
            if type(var) != dict:
                raise TypeError("interpolationDomain must be a dict containing"\
                                " the keys and the values of {'surfacePolygon','boundaryPolygon'}."\
                                "It is %s." % type(var))
            
            # check through data inside the dict
            for key in interpolationDomainKeys:
                
                if key not in var.keys():
                    raise TypeError("interpolationDomain must be a dict containing"\
                                    " the keys and the values of {'surfacePolygon','boundaryPolygon'}."\
                                    " It is %s." % type(var))
                
                # If none, don't check. 
                if var[key] is not None:
                    # Check data type
                    if type(var[key]) != _numpy.ndarray or var[key].dtype != float:
                        raise TypeError("interpolationDomain['%s'] must be float numpy array"\
                                        "It is %s(%s)" % (key,type(var[key]),var[key].dtype))
                    # Check data shape (:,2)
                    if var[key].shape[1] != 2:
                        raise ValueError("interpolationDomain[%s] must be of shape (nPoints,2)."\
                                         "It has shape %s." % (key, str(var[key].shape)))
               
            # Everything okay
            self.__interpolationDomain = var
                                
        #---------------------------------------------------------------------

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