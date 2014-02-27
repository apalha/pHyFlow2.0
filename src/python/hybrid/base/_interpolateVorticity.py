# -*- coding: utf-8 -*-
"""
Routines for coupling navier-stokes grid solver and the vortex particles.
Coupling done by transfering the grid vorticity to vortex particles.

:First Added:   2013-08-28
:Last Modified: 2013-10-04
:Copyright:     Copyright (C) 2013 Lento Manickathan, Artur Palha, **pHyFlow**
:Licence:       GNU GPL version 3 or any later version

"""
# Copyright (C) 2013 Lento Manickathan                                                               
#                                                                                                                                   
# This file is part of pHyFlow.                                                                                                      
#                                                                                                                                   
# pHyFlow is free software: you can redistribute it and/or modify                                                                    
# it under the terms of the GNU Lesser General Public License as published by                                                       
# the Free Software Foundation, either version 3 of the License, or                                                                 
# (at your option) any later version.                                                                                               
#                                                                                                                                   
# pHyFlow is distributed in the hope that it will be useful,                                                                         
# but WITHOUT ANY WARRANTY; without even the implied warranty of                                                                    
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the                                                                      
# GNU Lesser General Public License for more details.                                                                               
#                                                                                                                                   
# You should have received a copy of the GNU Lesser General Public License                                                          
# along with pHyFlow. If not, see <http://www.gnu.org/licenses/>.                                                                    


import numpy
from matplotlib.nxutils import points_inside_poly

__all__ = ['interpolateVorticity']


def interpolateVorticity(xBlob,yBlob,wBlob,sigma,overlap,xBoundsBlob,
                         yBoundsBlob,xGrid,yGrid,vortGrid,thetaGrid,hGrid,
                         xyPolygon):               
    """
    Interpolate the vorticity from Finite Element interpolation grid of the
    navier-stokes domain to the vortex blobs. The interpolation method will
    generate blobs inside the interpolation grid (according to the remeshing)
    strategy and will interpolate the local **circulation** onto the blobs. If 
    vortex blobs are inside already inside the interpolation domain, 
    their circulation is set to zero. Later, population control should be used
    to remove them.
    
    Strategy
    --------
        
        #. Generate temporary blobs, uniformly around (near and inside) the 
           interpolation grid.
           
        #. Remove the temporary blobs that are outside `boundary of interest` 
           (:math:`{{x,y}_{Grid}}_{boundary}`) of the interpolation grid using
           :py:func:`matplotlib.nxutil.points_in_polygon` point in polygon
           finder.
           
        #. Interpolate the vorticity from the interpolation grid, inside 
           `boundary of interest` onto the temporary blobs inside the
           interpolation grid using **linear interpolation** of the four grid
           verticies.
           
        #. Locate the origin blobs inside the `domain of interest` in the 
           interpolation grid and set their circulation to zero. First, locate
           inside the minimum bounding box of interpolation grid; second, use 
           :py:func:`matplotlib.nxutil.points_in_polygon` to find inside the 
           boundary of interest. This is the efficient approach for large
           number of blobs
           
        #. Finally, get new blobs by concatenating the temporary blobs and the
           origin blobs. The zero-circulation blobs should be removed later
           with population control

    Usage
    -----
    .. code-block:: python
    
        xBlobsNew,yBlobsNew,wBlobsNew = interpolateVorticity(xBlob,yBlob,wBlob,
                                            sigma,overlap,xBoundsBlob,
                                            yBoundsBlob,xGrid,yGrid,vortGrid,
                                            thetaGrid,hGrid,xyPolygon)


    Parameters 
    ---------- 
    xBlob : ndarray (float64), shape (nBlobs,)
            the :math:`x`-coordinates of the vortex blobs in `global coordinate 
            system`. 
        
    yBlob : ndarray (float64), shape (nBlobs,)
            the :math:`y`-coordinates of the vortex blobs in `global coordinate 
            system`. 
        
    wBlob : ndarray (float64), shape (nBlobs,)
            the **circulations** associated to each of the vortex blobs.
                
    sigma : float64
            the core size of all the vortex blobs

    overlap : float64
              the overlap of the blobs
              
    xBoundsBlob : ndarray (float64), shape (2,)
                  the :math:`x` bounds of the domain where the particles were
                  distributed originally
                      
    yBoundsBlob : ndarray (float64), shape (2,)
                  the :math:`y` bounds of the domain where the particles were
                  distributed originally
        
    
    xGrid : ndarray (float64), shape (nyGrid,nxGrid)
            the :math:`x` coordinates of the interpolation grid in `global
            coordinate system`.
            
    yGrid : ndarray (float64), shape (nyGrid,nxGrid)
            the :math:`y` coordinates of the interpolation grid in `global
            coordinate system`.                        
                
    vortGrid : ndarray (float64), shape (nyGrid,nxGrid)
                the **vorticity** at the interpolation grid coordinates.

    thetaGrid : float64, unit (radians)
                the **anti-clockwise** rotation angle of the interpolation grid
                w.r.t to the global coordinate system.

    hGrid : ndarray (float64), shape (2,)
            the local :math:`x,y` step size :math:`\Delta x, \Delta y` of the
            interpolation grid.               
                                       
    xyPolygon : ndarray (float64), shape (nVertices,nDim)
                a sequency of :math:`x,y`- coordinates of the vertices of the
                polygon describing the boundaries starting point. For multiple                      
                polygon lines, just concatentate the coordinates after closing
                the sub-polygon loops.
        
                    
    Returns
    -------
    xBlobNew : ndarray (float64), shape (nBlobs,)
               the :math:`x` coordinates of the new vortex blobs in `global 
               coordinate system`. 
        
    yBlobNew : ndarray (float64), shape (nBlobs,)
               the :math:`y` coordinates of the new vortex blobs in `global  
               coordinate system`. 
        
    wBlobNew : ndarray (float64), shape (nBlobs,)
               the **circulations** associated to each of the new vortex
               blobs.

    :First Added:   2013-08-28
    :Last Modified: 2013-10-04
    :Copyright:     Copyright (C) 2013 Lento Manickathan, Artur Palha, **pHyFlow**
    :Licence:       GNU GPL version 3 or any later version
    
    """   
    
    # New set of blobs
    xBlob = xBlob.copy()    
    yBlob = yBlob.copy()    
    wBlob = wBlob.copy()    
    
    # Define blob parameters
    hBlob = sigma * overlap
    
    #######################################################################
    #   DETERMINING THE LIMITS OF THE INTERPOLATION GRID
    
    # Bounds of the interpolation grid - ('in global')        
    xGridMin = xGrid.min()
    xGridMax = xGrid.max()
    yGridMin = yGrid.min()
    yGridMax = yGrid.max()
    
    #######################################################################

    ##########################################################################
    #   GENERATE BLOBS AROUND THE INTERPOLATION GRID - THE ONES TO BE USED FOR INTERPOLATION       
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

    # Generate 1D grid of the blobs  - nBlobsTemp evenly spaced
    xBlobTemp = numpy.arange(numpy.floor((xGridMin - xBoundsBlob[0])/hBlob)-1,
                             numpy.ceil( (xGridMax - xBoundsBlob[0])/hBlob)+1)*hBlob + xBoundsBlob[0] + 0.5*hBlob #nBlobsTemp evenly spaced in x

    yBlobTemp = numpy.arange(numpy.floor((yGridMin - yBoundsBlob[0])/hBlob)-1,
                             numpy.ceil( (yGridMax - yBoundsBlob[0])/hBlob)+1)*hBlob + yBoundsBlob[0] + 0.5*hBlob #nBlobsTemp evenly spaced in y

    # Generate 2D grid of blobs
    xBlobTemp,yBlobTemp = numpy.meshgrid(xBlobTemp,yBlobTemp)                             

    # Flatten the blobs
    xBlobTemp = xBlobTemp.flatten()
    yBlobTemp = yBlobTemp.flatten()

    # Find the blobs inside the interpolation grid
    maskBlobTemp = points_inside_poly(numpy.array([xBlobTemp, yBlobTemp]).T,
                                      xyPolygon)

    # Use the particles that are only inside the interpolation grid
    xBlobTemp = xBlobTemp[maskBlobTemp]
    yBlobTemp = yBlobTemp[maskBlobTemp]

    #######################################################################

    #######################################################################
    #   INTERPOLATE THE SOLUTION ONTO BLOBTEMP INSIDE INTERPOLATION GRID

    # Distance of the particles from the grid origin
    LxBlobTemp =   (xBlobTemp-xGrid[0,0])*numpy.cos(thetaGrid) \
                 + (yBlobTemp-yGrid[0,0])*numpy.sin(thetaGrid)
                 
    LyBlobTemp = - (xBlobTemp-xGrid[0,0])*numpy.sin(thetaGrid) \
                 + (yBlobTemp-yGrid[0,0])*numpy.cos(thetaGrid)

    # Indexing the particles w.r.t to the origin
    xLIndex = numpy.int64(numpy.floor(LxBlobTemp/hGrid[0]))
    yLIndex = numpy.int64(numpy.floor(LyBlobTemp/hGrid[1]))
    
    # Linear interpolation of vorticity
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
    #   Linear Interpolation:
    #   .. math::
    #       \\omega_{blob}\\left(x^{\prime},y^{\\prime}\\right)={\\textstyle \\sum_{i=1}^{4}\
    #               {\\omega_{grid}}_i\\cdot h_i\\left(x^{\\prime},y^{\\prime}\\right)},
    #
    #           where:
    #               left( x^\prime-\Delta x \right)\cdot\left( y^\prime-\Delta y \right)}{\Delta x \Delta y},
    #               h_2=-\frac{x^\prime\left( y^\prime-\Delta y \right)}{\Delta x \Delta y}
    #               h_3=\frac{x^\prime \cdot y^\prime}{\Delta x \Delta y}
    #               h_4=-\frac{y^\prime\left( x^\prime-\Delta x \right)}{\Delta x \Delta y}

    # Interpolating using the hat-function    
    
    # Distance of each blob from it's lower left grid box
    xPrime = LxBlobTemp - xLIndex*hGrid[0]
    yPrime = LyBlobTemp - yLIndex*hGrid[1]

    # Linear Interpolation weights
    h1 =    (xPrime - hGrid[0]) * (yPrime-hGrid[1]) / (hGrid[0]*hGrid[1])
    h2 =       - xPrime         * (yPrime-hGrid[1]) / (hGrid[0]*hGrid[1])
    h3 =         xPrime         *     yPrime        / (hGrid[0]*hGrid[1])
    h4 =       - yPrime         * (xPrime-hGrid[0]) / (hGrid[0]*hGrid[1])
        
    # Interpolating the (vorticity*area) = circulation from the grid the blobs    
    wBlobTemp = (h1*vortGrid[yLIndex,xLIndex]       +\
                 h2*vortGrid[yLIndex,xLIndex+1]     +\
                 h3*vortGrid[yLIndex+1,xLIndex+1]   +\
                 h4*vortGrid[yLIndex+1,xLIndex]) * (hGrid[0]*hGrid[1]) # Circulation = vort * area

    ##########################################################################

    
    ##########################################################################
    # ADD THE NEWLY GENERATED BLOB (FROM THE INTERPOLATION) TO THE EXISTING ONES 
    
    # Remove the old blobs inside the interpolation domain
    
    # Locate the index of blobs inside the bounding box square
    BlobInsideBBOXIndex = numpy.where((xBlob > xGridMin) & 
                                      (xBlob < xGridMax) & 
                                      (yBlob > yGridMin) & 
                                      (yBlob < yGridMax))[0]
    

    maskBlob = points_inside_poly(numpy.array([xBlob[BlobInsideBBOXIndex],yBlob[BlobInsideBBOXIndex]]).T,
                                  xyPolygon)         
      

    """                                                 
    # Ensuring that the total change in circulation is zero
    totalCirculationOld = wBlob[BlobInsideBBOXIndex[maskBlob]].sum() # numpy.abs(wBlob[maskBlob]).sum()    
    print "total Circulation old : " + str(totalCirculationOld)
    
    totalCirculationNew = wBlobTemp.sum()                            # numpy.abs(wBlobTemp).sum() 
    print "total Circulation new : " + str(totalCirculationNew)
    
    dtotalCirculation   = totalCirculationOld - totalCirculationNew
    print "dGamma total : " + str(dtotalCirculation)
    
    #print "dGamma per blob : " + str(dtotalCirculation/wBlobTemp.shape[0])

    wBlobTempBefore = wBlobTemp.copy()
    #print "wBlobTemp before : " + str(wBlobTempBefore)
    
    #wBlobTemp += dtotalCirculation/wBlobTemp.shape[0]                    
    #print "wBlobTemp after : " + str(wBlobTempBefore)

    #wBlobTemp += dtotalCirculation * numpy.abs(wBlobTemp) / numpy.abs(wBlobTemp).sum()    
    
    print "wBlobTemp diff : " + str (wBlobTemp.sum() - wBlobTempBefore.sum()) + ", should be equal to dGamma total"
                    
    print "Ensure correction is working: dGamma = " + str(wBlobTemp.sum() - wBlob[BlobInsideBBOXIndex[maskBlob]].sum()) \
        + ", should be equal to zero."        
            
    """
           
    # Zero the blobs inside the interpolation grid
    wBlob[BlobInsideBBOXIndex[maskBlob]] = 0.0
    
    
    # Concatenated the blobs from vortex method (which are not inside) and the
    # newly generated blobs        
    xBlobNew = numpy.concatenate((xBlob,xBlobTemp))
    yBlobNew = numpy.concatenate((yBlob,yBlobTemp))
    wBlobNew = numpy.concatenate((wBlob,wBlobTemp))

    #wBlobNew += dtotalCirculation * numpy.abs(wBlobNew) / numpy.abs(wBlobNew).sum()
    
    #    print "total circulation: " + str(wBlobNew.sum())
    #    xBlobNew = xBlob.copy()
    #    yBlobNew = yBlob.copy()
    #    wBlobNew = wBlob.copy()
        
    # Return the (new+old) vortex blobs
    return xBlobNew,yBlobNew,wBlobNew
