"""

  General routines for plotting vortex blob related quantities: velocity field
  vorticity field, etc.

"""
# Copyright (C) 2013 Artur Palha                                                                                                     
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
#                                                                                                                                   
# First added:  2013-05-22                                                                                                          
# Last changed: 2013-05-28
# -*- coding: utf-8 -*-

__all__ = ['plotvorticity','plotvelocity','plotblobs']

import numpy
import pylab
import dolfin
from pHyFlow import options

from pHyFlow.vortex.base.induced import velocity
from pHyFlow.vortex.base.induced import vorticity,vorticity_blobs

# set interactive mode on
pylab.ion()


#def plot(xBounds,yBounds,xBlob,yBlob,wBlob,sigma,overlap=0.5,k=2,kernel=1,\
#         nPlotPoints=[20,20],hardware=0,blocksize=128,method=1,field=0,wType=0,figureHandle=None):
#    """
#        Plot vortex blobs induced fields (velocity and vorticity) in a region.
#        
#        Usage
#        -----
#            figureHandle = plot(xBounds,yBounds,xBlob,yBlob,wBlob,
#                                sigma,k=2,kernel=1,nPlotPoints=20,
#                                hardware=0,blocksize=128,method=1,
#                                field=0,wType=0,figureHandle=None)
#            
#        Parameters
#        ----------
#            xBounds :: the x bounds of the box where to plot the induced field
#            -------    (type: numpy.ndarray (float64); shape: (2,))
#            
#            yBounds :: the y bounds of the box where to plot the induced field
#            -------    (type: numpy.ndarray (float64); shape: (2,))
#            
#            xBlob :: the x coordinates of the vortex blobs
#            -----    (type: numpy.ndarray (float64); shape: (nBlobs,))
#                    
#            yBlob :: the y coordinates of the vortex blobs
#            -----    (type: numpy.ndarray (float64); shape: (nBlobs,))
#                     
#            wBlob :: the circulations associated to each of the vortex blobs
#            -----    (type: numpy.ndarray (float64); shape: (nBlobs,))
#                     
#            sigma :: the core size of all the vortex blobs
#            -----    (type: float (64bits); shape: single value)
#            
#            overlap :: the overlap of the blobs
#                       (default value is 0.5)
#                       (type: float (64bits); shape: single value)
#                     
#            k :: the core size multiplication constant of all the vortex blobs
#                 typical values are 1,2,4 (default value is 2)
#                 (type: float (64bits); shape: single value)
#                 
#            kernel :: the type of kernel of all the vortex blobs
#                      available kernels are: 0 ('cutoff'), 1 ('gauss')
#                      (default value is gauss kernel)
#                      (type: int; shape: single value)
#                      
#            nPlotPoints :: the number of points in x and y directions to use
#                           to generate the plotting grid
#                           (default value is [20,20])
#                           (type: numpy.ndarray (int32); shape: (2,))
#                                
#            hardware :: the hardware to use to compute the induced fields
#                        can be the 0 (for CPU) or 1 (for GPU)
#                        (default value is CPU)
#                        (type: int; shape: single value)
#                        
#            blocksize :: the size of gpu memory block size
#                         (default value 128)
#                         (type: int; shape: single value)
#            
#            method :: the method used to compute the induced fields can be
#                      0 (for FMM) or 1 (for direct calculation)
#                      (default value is direct calculation)
#                      (type: int; shape: single value)
#                      
#            field :: the induced field to plot can be 0 (vorticity) or
#                     1 (velocity)
#                     (default value is 0)
#                     (type: int; shape: single value)
#                     
#            wType :: the method used to compute the vorticity field
#                     0 (vorticity obtained from Gamma*/(deltaX*deltaY)) or
#                     1 (vorticity obtained by evaluating the vorticity functions)
#                     option 0 is faster but less detailed as the spacing of the
#                     plot points is equal to the cell size associated to each
#                     blob
#                     (default value is 0)
#                     (type: int; shape: single value)
#                     
#            figureHandle :: the handle for the figure of the plot, if no value
#                            is given a new figure is created, if a value is given
#                            use that figure to plot
#                            (default value is None)
#                            (type: int; shape: single value)
#          
#                     
#        Returns
#        -------
#            figureHandle :: the handle for the figure of the plot
#                            (type: int; shape: single value)
#                  
#                  
#        First added:     2013-06-27
#
#        Copyright (C) 2013 Artur Palha
#                           pHyFlow
#    """
#    
#    """    
#        Reviews:   
#    """
#    
#    if field == options.W_FIELD: # plot vorticity field
#        figureHandle = _plotvorticity(xBounds,yBounds,\
#                                      xBlob,yBlob,wBlob,sigma,overlap,k,kernel,\
#                                      nPlotPoints,hardware,blocksize,method,wType,figureHandle)
#    elif field == options.V_FIELD: # plot velocity field
#        print 'VELOCITY FIELD'        
#        figureHandle = 1
#        
#    # return figure handle to the plot
#    return figureHandle
#    

def plotblobs(xBoundsPlot,yBoundsPlot,xBlob,yBlob,wBlob,figureHandle=None,title=None):
    """
        Plot vortex blobs in a region.
        
        Usage
        -----
            figureHandle = plotblobs(xBoundsPlot,yBoundsPlot,xBlob,yBlob,wBlob,
                                     figureHandle=None,title=None)
            
        Parameters
        ----------
            xBoundsPlot :: the x bounds of the box where to plot the induced field
            -----------   (type: numpy.ndarray (float64); shape: (2,))
            
            yBoundsPlot :: the y bounds of the box where to plot the induced field
            -----------    (type: numpy.ndarray (float64); shape: (2,))
            
            xBlob :: the x coordinates of the vortex blobs
            -----    (type: numpy.ndarray (float64); shape: (nBlobs,))
                    
            yBlob :: the y coordinates of the vortex blobs
            -----    (type: numpy.ndarray (float64); shape: (nBlobs,))
                     
            wBlob :: the circulations associated to each of the vortex blobs
            -----    (type: numpy.ndarray (float64); shape: (nBlobs,))
                     
            figureHandle :: the handle for the figure of the plot, if no value
                            is given a new figure is created, if a value is given
                            use that figure to plot
                            (default value is None)
                            (type: int; shape: single value)
                            
            title :: the title of the plot
                     (default values is None)
                     (type: str; shape: the shape of a string (nchar,))                            
                     
        Returns
        -------
            figureHandle :: the handle for the figure of the plot
                            (type: int; shape: single value)
                  
                  
        First added:     2013-07-22

        Copyright (C) 2013 Artur Palha
                           pHyFlow
    """
    
    """    
        Reviews:   
    """
    # generate a figure
    if figureHandle == None: # if no figureHandle is given make a new figure
        figureHandle = pylab.figure()
    else: # if a figure handle is given, use that figure to plot
        figureHandle = pylab.figure(figureHandle)
        # clear the figure    
        pylab.clf()
    
    # make the scatter plot
    pylab.scatter(xBlob,yBlob,c=wBlob)
    
    # set the axis and title and draw
    pylab.axis('scaled')    
    pylab.axis(numpy.concatenate((xBoundsPlot,yBoundsPlot)))
    pylab.title(title)                         
    pylab.draw()

    return figureHandle
    
    
    
def plotvorticity(xBounds,yBounds,xBlob,yBlob,wBlob,sigma,overlap,k=2,kernel=1,\
                   nPlotPoints=[20,20],hardware=0,blocksize=128,method=1,wType=0,figureHandle=None,title=None,plotType=0,interactive=False):
    """
        Plot vortex blobs induced fields (velocity and vorticity) in a region.
        
        Usage
        -----
            figureHandle = plotvorticity(xBounds,yBounds,xBlob,yBlob,wBlob,
                                         sigma,k=2,kernel=1,nPlotPoints=[20,20],
                                         hardware=0,blocksize=128,method=1,
                                         wType=0,figureHandle=None,title=None,
                                         plotType=0,interactive=False)
            
        Parameters
        ----------
            xBounds :: the x bounds of the box where to plot the induced field
            -------    (type: numpy.ndarray (float64); shape: (2,))
            
            yBounds :: the y bounds of the box where to plot the induced field
            -------    (type: numpy.ndarray (float64); shape: (2,))
            
            xBlob :: the x coordinates of the vortex blobs
            -----    (type: numpy.ndarray (float64); shape: (nBlobs,))
                    
            yBlob :: the y coordinates of the vortex blobs
            -----    (type: numpy.ndarray (float64); shape: (nBlobs,))
                     
            wBlob :: the circulations associated to each of the vortex blobs
            -----    (type: numpy.ndarray (float64); shape: (nBlobs,))
                     
            sigma :: the core size of all the vortex blobs
            -----    (type: float (64bits); shape: single value)
            
            overlap :: the overlap of the blobs
                       (default value is 0.5)
                       (type: float (64bits); shape: single value)
                     
            k :: the core size multiplication constant of all the vortex blobs
                 typical values are 1,2,4 (default value is 2)
                 (type: float (64bits); shape: single value)
                 
            kernel :: the type of kernel of all the vortex blobs
                      available kernels are: 0 ('cutoff'), 1 ('gauss')
                      (default value is gauss kernel)
                      (type: int; shape: single value)
                      
            nPlotPoints :: the number of points in x and y directions to use
                           to generate the plotting grid
                           (default value is [20,20])
                           (type: numpy.ndarray (int32); shape: (2,))
                                
            hardware :: the hardware to use to compute the induced fields
                        can be the 0 (for CPU) or 1 (for GPU)
                        (default value is CPU)
                        (type: int; shape: single value)
                        
            blocksize :: the size of gpu memory block size
                         (default value 128)
                         (type: int; shape: single value)
            
            method :: the method used to compute the induced fields can be
                      0 (for FMM) or 1 (for direct calculation)
                      (default value is direct calculation)
                      (type: int; shape: single value)
                     
            wType :: the method used to compute the vorticity field
                     0 (vorticity obtained from Gamma*/(deltaX*deltaY)) or
                     1 (vorticity obtained by evaluating the vorticity functions)
                     option 0 is faster but less detailed as the spacing of the
                     plot points is equal to the cell size associated to each
                     blob
                     (default value is 0)
                     (type: int; shape: single value)
                     
            figureHandle :: the handle for the figure of the plot, if no value
                            is given a new figure is created, if a value is given
                            use that figure to plot
                            (default value is None)
                            (type: int; shape: single value)
                            
            title :: the title of the plot
                     (default values is None)
                     (type: str; shape: the shape of a string (nchar,))                            
                            
            plotType :: define the type of plot
                        0 - pylab plot
                        1 - dolfin plot
                        (default value is 0)
                        (type: int; shape: single value)                            
            
            interactive :: define if plot is interactive (True) or not (False)
                           only used for plotting with dolfin
                           (default value is False)
                           (type: bool; shape: single value)
          
                     
        Returns
        -------
            figureHandle :: the handle for the figure of the plot
                            (type: int; shape: single value)
                  
                  
        First added:     2013-06-27

        Copyright (C) 2013 Artur Palha
                           pHyFlow
    """
    
    """    
        Reviews:   1- Changed to a public function. (2013-07-04)
                   1- Added title option. (2013-07-16)
    """    
    
    if wType == options.WBLOB_VORTICITY: # use blobs circulation
        if plotType == options.PYLAB_PLOT: # use pylab plot
            figureHandle = _plotvorticity_blobs_pylab(xBounds,yBounds,xBlob,yBlob,wBlob,sigma,overlap,figureHandle,title)
        
        if plotType == options.DOLFIN_PLOT: # use dolfin plot
            figureHandle = _plotvorticity_blobs_dolfin(xBounds,yBounds,xBlob,yBlob,wBlob,sigma,overlap,figureHandle,title,interactive)
        
    elif wType == options.WFUNCTION_VORTICITY: # use blobs function
        if plotType == options.PYLAB_PLOT: # use pylab plot                
            figureHandle = _plotvorticity_function_pylab(xBounds,yBounds,\
                                                         xBlob,yBlob,wBlob,sigma,k,kernel,\
                                                         nPlotPoints,hardware,blocksize,method,figureHandle,title)
                                                         
        if plotType == options.DOLFIN_PLOT: # use dolfin plot                
            figureHandle = _plotvorticity_function_dolfin(xBounds,yBounds,\
                                                          xBlob,yBlob,wBlob,sigma,k,kernel,\
                                                          nPlotPoints,hardware,blocksize,method,figureHandle,title,interactive)                                                         
    
    # return figure handle with the plot                                               
    return figureHandle
    
    
    
def _plotvorticity_blobs_pylab(xBounds,yBounds,xBlob,yBlob,wBlob,sigma,overlap,figureHandle,title):
    """
        Plot (with pylab) vortex blobs induced vorticity in a region using blobs
        circulation.
        
        Usage
        -----
            figureHandle = _plotvorticity_blobs_pylab(xBounds,yBounds,
                                                xBlob,yBlob,wBlob,
                                                sigma,overlap, figureHandle,title)
            
        Parameters
        ----------
            xBounds :: the x bounds of the box where to plot the induced field
            -------    (type: numpy.ndarray (float64); shape: (2,))
            
            yBounds :: the y bounds of the box where to plot the induced field
            -------    (type: numpy.ndarray (float64); shape: (2,))
            
            xBlob :: the x coordinates of the vortex blobs
            -----    (type: numpy.ndarray (float64); shape: (nBlobs,))
                    
            yBlob :: the y coordinates of the vortex blobs
            -----    (type: numpy.ndarray (float64); shape: (nBlobs,))
                     
            wBlob :: the circulations associated to each of the vortex blobs
            -----    (type: numpy.ndarray (float64); shape: (nBlobs,))
                     
            sigma :: the core size of all the vortex blobs
            -----    (type: float (64bits); shape: single value)
            
            overlap :: the overlap of the blobs
            -------    (type: float (64bits); shape: single value)
                       
           figureHandle :: the number of the figure where to plot. If None, a
           ------------    a new figure is created.
                           (type: int; shape: single value)
                           
           title :: the title of the plot
           -----    (default values is None)
                    (type: str; shape: the shape of a string (nchar,)) 

                     
        Returns
        -------
            figureHandle :: the handle for the figure that holds the plot
                            (type: int; shape: single value)
                  
                  
        First added:     2013-06-27

        Copyright (C) 2013 Artur Palha
                           pHyFlow
    """
    
    """    
        Reviews:   1- Changed name from _plotvorticity_blob to 
                      _plotvorticity_blobs_pylab in order to make another
                      function where plotting is done with dolfin. (2013-07-04)
    """    
    
    # compute the vorticity associated to each blob, as given by:
    #   w = Gamma * deltaX * deltaY
    # compute also the coordinates of the blobs
    xW,yW,w = vorticity_blobs(xBlob,yBlob,wBlob,sigma,overlap,xBounds,yBounds)

    # generate a figure
    if figureHandle == None: # if no figureHandle is given make a new figure
        plotFigure = pylab.figure()
    else: # if a figure handle is given, use that figure to plot
        plotFigure = pylab.figure(figureHandle)
    
    # plot the blobs
    pylab.imshow(w,extent=(xW.min()-0.5*sigma*overlap,\
                           xW.max()+0.5*sigma*overlap,\
                           yW.min()-0.5*sigma*overlap,\
                           yW.max()+0.5*sigma*overlap),\
                           interpolation='nearest',origin='lower')
    
    # set the title of the plot
    pylab.title(title)
    
    # draw the figure
    pylab.draw()
    
    # return figure handle to the plot
    return plotFigure.number



def _plotvorticity_blobs_dolfin(xBounds,yBounds,xBlob,yBlob,wBlob,sigma,overlap,figureHandle,title,interactive):
    """
        Plot (with dolfin) vortex blobs induced vorticity in a region using blobs
        circulation.
        
        Usage
        -----
            figureHandle = _plotvorticity_blobs_dolfin(xBounds,yBounds,
                                                xBlob,yBlob,wBlob,
                                                sigma,overlap, figureHandle,title,
                                                interactive)
            
        Parameters
        ----------
            xBounds :: the x bounds of the box where to plot the induced field
            -------    (type: numpy.ndarray (float64); shape: (2,))
            
            yBounds :: the y bounds of the box where to plot the induced field
            -------    (type: numpy.ndarray (float64); shape: (2,))
            
            xBlob :: the x coordinates of the vortex blobs
            -----    (type: numpy.ndarray (float64); shape: (nBlobs,))
                    
            yBlob :: the y coordinates of the vortex blobs
            -----    (type: numpy.ndarray (float64); shape: (nBlobs,))
                     
            wBlob :: the circulations associated to each of the vortex blobs
            -----    (type: numpy.ndarray (float64); shape: (nBlobs,))
                     
            sigma :: the core size of all the vortex blobs
            -----    (type: float (64bits); shape: single value)
            
            overlap :: the overlap of the blobs
            -------    (type: float (64bits); shape: single value)
                       
           figureHandle :: the handle for the figure that holds the plot
           ------------    in this case a dolfin plot
                           (type: dolfin.cpp.io.VTKPlotter; shape: single value)
                           
           title :: the title of the plot
           -----    (default values is None)
                    (type: str; shape: the shape of a string (nchar,))                           
                           
           interactive :: define if plot is interactive (True) or not (False)
           -----------    (type: bool; shape: single value)

                     
        Returns
        -------
            figureHandle :: the handle for the figure that holds the plot
                            in this case a dolfin plot
                            (type: dolfin.cpp.io.VTKPlotter; shape: single value)
                  
                  
        First added:     2013-07-04

        Copyright (C) 2013 Artur Palha
                           pHyFlow
    """
    
    """    
        Reviews:   
    """    
    
    # compute the vorticity associated to each blob, as given by:
    #   w = Gamma * deltaX * deltaY
    # compute also the coordinates of the blobs
    xW,yW,w = vorticity_blobs(xBlob,yBlob,wBlob,sigma,overlap,xBounds,yBounds)

    # generate the mesh where to plot the data
    # this mesh extends from (xWmin,yWmin) to (xWmax,yWmax) and has nxBlobs points
    # in the x direction and nyBlobs points in the y direction
    mesh = dolfin.RectangleMesh(xW.min(),yW.min(),xW.max(),yW.max(),xW.shape[1]-1,xW.shape[0]-1)
    
    # generate the function space: linear Lagrange interpolants
    V = dolfin.FunctionSpace(mesh,'CG',1)
    
    # get the matrix that converts degrees of freedom from the mesh into
    # finite element degrees of freedom
    dof_to_vertex = V.dofmap().dof_to_vertex_map(mesh)
    
    # generate the function that will store the vorticity to plot
    wPlot = dolfin.Function(V)
    
    # store the vorticity in this function
    wPlot.vector()[dof_to_vertex] = w.flatten()
    
    # generate the plot
    if figureHandle == None: # if no figureHandle is given make a new figure
        figureHandle = dolfin.plot(wPlot,title=title)
    else: # if a figure handle is given, use that figure to plot
        figureHandle.plot(wPlot)
        
    # return figure handle to the plot
    return figureHandle
    
    
    
def _plotvorticity_function_pylab(xBounds,yBounds,xBlob,yBlob,wBlob,sigma,k,kernel,\
                            nPlotPoints,hardware,blocksize,method,figureHandle,title):
    """
        Plot vortex blobs induced vorticity in a region using vorticity function.
        
        Usage
        -----
            figureHandle = _plotvorticity_function_pylab(xBounds,yBounds,
                                                   xBlob,yBlob,wBlob,sigma,k,kernel,\
                                                   nPlotPoints,hardware,blocksize,\
                                                   method, figureHandle,title)
            
        Parameters
        ----------
            xBounds :: the x bounds of the box where to plot the induced field
            -------    (type: numpy.ndarray (float64); shape: (2,))
            
            yBounds :: the y bounds of the box where to plot the induced field
            -------    (type: numpy.ndarray (float64); shape: (2,))
            
            xBlob :: the x coordinates of the vortex blobs
            -----    (type: numpy.ndarray (float64); shape: (nBlobs,))
                    
            yBlob :: the y coordinates of the vortex blobs
            -----    (type: numpy.ndarray (float64); shape: (nBlobs,))
                     
            wBlob :: the circulations associated to each of the vortex blobs
            -----    (type: numpy.ndarray (float64); shape: (nBlobs,))
                     
            sigma :: the core size of all the vortex blobs
            -----    (type: float (64bits); shape: single value)
            
            k :: the core size multiplication constant of all the vortex blobs
            -    typical values are 1,2,4 (default value is 2)
                 (type: float (64bits); shape: single value)
                 
            kernel :: the type of kernel of all the vortex blobs
            ------    available kernels are: 0 ('cutoff'), 1 ('gauss')
                      (default value is gauss kernel)
                      (type: int; shape: single value)
                      
            nPlotPoints :: the number of points in x and y directions to use
            -----------    to generate the plotting grid
                           (default value is [20,20])
                           (type: numpy.ndarray (int32); shape: (2,))
                                
            hardware :: the hardware to use to compute the induced fields
            --------    can be the 0 (for CPU) or 1 (for GPU)
                        (default value is CPU)
                        (type: int; shape: single value)
                        
            blocksize :: the size of gpu memory block size
            ---------    (default value 128)
                         (type: int; shape: single value)
            
            method :: the method used to compute the induced fields can be
            ------    0 (for FMM) or 1 (for direct calculation)
                      (default value is direct calculation)
                      (type: int; shape: single value)
                      
            figureHandle :: the number of the figure where to plot. If None, a
            ------------    a new figure is created.
                           (type: int; shape: single value)

            title :: the title of the plot
            -----    (default values is None)
                    (type: str; shape: the shape of a string (nchar,)) 

                     
        Returns
        -------
            figureHandle :: the handle for the figure that holds the plot
                            (type: int; shape: single value)
                  
                  
        First added:     2013-07-08

        Copyright (C) 2013 Artur Palha
                           pHyFlow
    """
    
    """    
        Reviews:   
    """    

    # generate the points where to plot vorticity
    # generate evenly spaced points in a rectangular box of dimensions
    # (xBounds[1]-xBounds[0]) in x direction and (yBounds[1]-yBounds[0]) in y direction
    # with nPlotPoints[0] points in x direction and nPlotPoints[1] points in y direction
    xW,yW = numpy.meshgrid(numpy.linspace(xBounds[0],xBounds[1],nPlotPoints[0]),\
                           numpy.linspace(yBounds[0],yBounds[1],nPlotPoints[1]))

    # compute the x and y spacing of plotting points
    deltaXplot = (xBounds[1]-xBounds[0])/(nPlotPoints[0]-1)
    deltaYplot = (yBounds[1]-yBounds[0])/(nPlotPoints[1]-1)
    
    # compute the induced vorticity at each of the plot points
    w = vorticity(xBlob,yBlob,wBlob,sigma,k=2,kernel=1,\
                  xEval=xW.flatten(),yEval=yW.flatten(),hardware=hardware,\
                  blocksize=blocksize,method=method)

    # generate a figure
    if figureHandle == None: # if no figureHandle is given make a new figure
        plotFigure = pylab.figure()
    else: # if a figure handle is given, use that figure to plot
        plotFigure = pylab.figure(figureHandle)
    
    # plot the blobs
    pylab.imshow(w.reshape([nPlotPoints[1],nPlotPoints[0]]),extent=(xW.min()-0.5*deltaXplot,\
                           xW.max()+0.5*deltaXplot,\
                           yW.min()-0.5*deltaYplot,\
                           yW.max()+0.5*deltaYplot),\
                           interpolation='nearest',origin='lower')

    # add the tile
    pylab.title(title)    
    
    # draw the figure
    pylab.draw()
    
    # return figure handle to the plot
    return plotFigure.number
    
    
    
def _plotvorticity_function_dolfin(xBounds,yBounds,xBlob,yBlob,wBlob,sigma,k,kernel,\
                            nPlotPoints,hardware,blocksize,method,figureHandle,title,interactive):
    """
        Plot vortex blobs induced vorticity in a region using vorticity function.
        
        Usage
        -----
            figureHandle = _plotvorticity_function_dolfin(xBounds,yBounds,
                                                   xBlob,yBlob,wBlob,sigma,k,kernel,\
                                                   nPlotPoints,hardware,blocksize,method,\
                                                   figureHandle,title,interactive)
            
        Parameters
        ----------
            xBounds :: the x bounds of the box where to plot the induced field
            -------    (type: numpy.ndarray (float64); shape: (2,))
            
            yBounds :: the y bounds of the box where to plot the induced field
            -------    (type: numpy.ndarray (float64); shape: (2,))
            
            xBlob :: the x coordinates of the vortex blobs
            -----    (type: numpy.ndarray (float64); shape: (nBlobs,))
                    
            yBlob :: the y coordinates of the vortex blobs
            -----    (type: numpy.ndarray (float64); shape: (nBlobs,))
                     
            wBlob :: the circulations associated to each of the vortex blobs
            -----    (type: numpy.ndarray (float64); shape: (nBlobs,))
                     
            sigma :: the core size of all the vortex blobs
            -----    (type: float (64bits); shape: single value)
            
            k :: the core size multiplication constant of all the vortex blobs
            -    typical values are 1,2,4 (default value is 2)
                 (type: float (64bits); shape: single value)
                 
            kernel :: the type of kernel of all the vortex blobs
            ------    available kernels are: 0 ('cutoff'), 1 ('gauss')
                      (default value is gauss kernel)
                      (type: int; shape: single value)
                      
            nPlotPoints :: the number of points in x and y directions to use
            -----------    to generate the plotting grid
                           (default value is [20,20])
                           (type: numpy.ndarray (int32); shape: (2,))
                                
            hardware :: the hardware to use to compute the induced fields
            --------    can be the 0 (for CPU) or 1 (for GPU)
                        (default value is CPU)
                        (type: int; shape: single value)
                        
            blocksize :: the size of gpu memory block size
            ---------    (default value 128)
                         (type: int; shape: single value)
            
            method :: the method used to compute the induced fields can be
            ------    0 (for FMM) or 1 (for direct calculation)
                      (default value is direct calculation)
                      (type: int; shape: single value)
                      
            figureHandle :: the handle for the figure that holds the plot
            ------------    in this case a dolfin plot
                            (type: dolfin.cpp.io.VTKPlotter; shape: single value)
                            
            title :: the title of the plot
            -----    (default values is None)
                    (type: str; shape: the shape of a string (nchar,)) 
                            
            interactive :: define if plot is interactive (True) or not (False)
            -----------    (type: bool; shape: single value)
                     
        Returns
        -------
            figureHandle :: the handle for the figure that holds the plot
                            in this case a dolfin plot
                            (type: dolfin.cpp.io.VTKPlotter; shape: single value)
                  
                  
        First added:     2013-07-04

        Copyright (C) 2013 Artur Palha
                           pHyFlow
    """
    
    """    
        Reviews:   
    """    
    
    # generate the points where to plot vorticity
    # generate evenly spaced points in a rectangular box of dimensions
    # (xBounds[1]-xBounds[0]) in x direction and (yBounds[1]-yBounds[0]) in y direction
    # with nPlotPoints[0] points in x direction and nPlotPoints[1] points in y direction
    xW,yW = numpy.meshgrid(numpy.linspace(xBounds[0],xBounds[1],nPlotPoints[0]),\
                           numpy.linspace(yBounds[0],yBounds[1],nPlotPoints[1]))

    # compute the induced vorticity at each of the plot points
    w = vorticity(xBlob,yBlob,wBlob,sigma,k=2,kernel=1,\
                  xEval=xW.flatten(),yEval=yW.flatten(),hardware=hardware,\
                  blocksize=blocksize,method=method)
                  
    # generate the mesh where to plot the data
    # this mesh extends from (xWmin,yWmin) to (xWmax,yWmax) and has nxBlobs points
    # in the x direction and nyBlobs points in the y direction
    mesh = dolfin.RectangleMesh(xW.min(),yW.min(),xW.max(),yW.max(),xW.shape[1]-1,xW.shape[0]-1)
    
    # generate the function space: linear Lagrange interpolants
    V = dolfin.FunctionSpace(mesh,'CG',1)
    
    # get the matrix that converts degrees of freedom from the mesh into
    # finite element degrees of freedom
    dof_to_vertex = V.dofmap().dof_to_vertex_map(mesh)
    
    # generate the function that will store the vorticity to plot
    wPlot = dolfin.Function(V)
    
    # store the vorticity in this function
    wPlot.vector()[dof_to_vertex] = w
    
    # generate the plot
    if figureHandle == None: # if no figureHandle is given make a new figure
        figureHandle = dolfin.plot(wPlot,title=title)
    else: # if a figure handle is given, use that figure to plot
        figureHandle.plot(wPlot)
        
    # return figure handle to the plot
    return figureHandle
    
    

def plotvelocity(xBounds,yBounds,xBlob,yBlob,wBlob,sigma,overlap,k=2,kernel=1,\
                   nPlotPoints=20,hardware=0,blocksize=128,method=1,figureHandle=None):
    """
        Plot vortex blobs induced fields (velocity and vorticity) in a region.
        
        Usage
        -----
            figureHandle = plotvorticity(xBounds,yBounds,xBlob,yBlob,wBlob,
                                         sigma,k=2,kernel=1,nPlotPoints=20,
                                         hardware=0,blocksize=128,method=1,
                                         wType=0,figureHandle=None)
            
        Parameters
        ----------
            xBounds :: the x bounds of the box where to plot the induced field
            -------    (type: numpy.ndarray (float64); shape: (2,))
            
            yBounds :: the y bounds of the box where to plot the induced field
            -------    (type: numpy.ndarray (float64); shape: (2,))
            
            xBlob :: the x coordinates of the vortex blobs
            -----    (type: numpy.ndarray (float64); shape: (nBlobs,))
                    
            yBlob :: the y coordinates of the vortex blobs
            -----    (type: numpy.ndarray (float64); shape: (nBlobs,))
                     
            wBlob :: the circulations associated to each of the vortex blobs
            -----    (type: numpy.ndarray (float64); shape: (nBlobs,))
                     
            sigma :: the core size of all the vortex blobs
            -----    (type: float (64bits); shape: single value)
            
            overlap :: the overlap of the blobs
                       (default value is 0.5)
                       (type: float (64bits); shape: single value)
                     
            k :: the core size multiplication constant of all the vortex blobs
                 typical values are 1,2,4 (default value is 2)
                 (type: float (64bits); shape: single value)
                 
            kernel :: the type of kernel of all the vortex blobs
                      available kernels are: 0 ('cutoff'), 1 ('gauss')
                      (default value is gauss kernel)
                      (type: int; shape: single value)
                      
            nPlotPoints :: the number of points in x and y directions to use
                           to generate the plotting grid
                           (default value is [20,20])
                           (type: numpy.ndarray (int32); shape: (2,))
                                
            hardware :: the hardware to use to compute the induced fields
                        can be the 0 (for CPU) or 1 (for GPU)
                        (default value is CPU)
                        (type: int; shape: single value)
                        
            blocksize :: the size of gpu memory block size
                         (default value 128)
                         (type: int; shape: single value)
            
            method :: the method used to compute the induced fields can be
                      0 (for FMM) or 1 (for direct calculation)
                      (default value is direct calculation)
                      (type: int; shape: single value)
                     
            figureHandle :: the handle for the figure of the plot, if no value
                            is given a new figure is created, if a value is given
                            use that figure to plot
                            (default value is None)
                            (type: int; shape: single value)
          
                     
        Returns
        -------
            figureHandle :: the handle for the figure of the plot
                            (type: int; shape: single value)
                  
                  
        First added:     2013-07-04

        Copyright (C) 2013 Artur Palha
                           pHyFlow
    """
    
    """    
        Reviews:
    """    
    
    print 'VELOCITY PLOT'
    
    # return figure handle with the plot                                               
    return 1
    
    
