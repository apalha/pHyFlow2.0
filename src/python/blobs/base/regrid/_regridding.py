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
# First added:  2013-07-10                                                                                                          
# Last changed: 2013-07-15
# -*- coding: utf-8 -*-

__all__ = ['Regrid','PopulationControl']

import numpy
from pHyFlow.blobs import blobOptions
from pHyFlow.cpp.blobs.remesh import *

import scipy.sparse 

def Regrid(xBlob,yBlob,wBlob,sigma,overlap,xBounds,yBounds,interpKernel=0,c=0.0):
    """
        Redistribute vortex blobs over evenly spaced grid points with separation
        h = sigma * overlap.
        The interpolation kernel can be one of the following:
           - M4' (0)
        
        Usage
        -----
            xBlobNew,yBlobNew,wBlobNew = regrid(xBlob,yBlob,wBlob,sigma,overlap
                                                xBounds,yBounds,interpKernel=0)
            
        Parameters
        ----------
            xBlob :: the x coordinates of the vortex blobs
            -----    (type: numpy.ndarray (float64); shape: (nBlobs,))
                    
            yBlob :: the y coordinates of the vortex blobs
            -----    (type: numpy.ndarray (float64); shape: (nBlobs,))
                     
            wBlob :: the circulations associated to each of the vortex blobs
            -----    (type: numpy.ndarray (float64); shape: (nBlobs,))
                     
            sigma :: the core size of all the vortex blobs
            -----    (type: float (64bits); shape: single value)
            
            overlap :: the overlap of the blobs
            -------   (default value is 0.5)
                      (type: float (64bits); shape: single value)
                       
            xBounds :: the x bounds of the domain where the particles were
            -------    distributed originally
                       (type: float (64bits); shape: (2,))
            
            yBounds :: the y bounds of the domain where the particles were
            -------    distributed originally
                       (type: float (64bits); shape: (2,))
                       
            interpKernel :: the interpolating kernel using for remeshing
                            (default value 0 (M4'))
                            (type: int (64bits); shape: single value)
          
                     
        Returns
        -------
            xBlobNew :: the x coordinates of the new vortex blobs
            --------    (type: numpy.ndarray (float64); shape: (nBlobsNew,))
                    
            yBlobNew :: the y coordinates of the new vortex blobs
            --------    (type: numpy.ndarray (float64); shape: (nBlobsNew,))
                     
            wBlobNew :: the circulations associated to each of the new vortex blobs
            --------    (type: numpy.ndarray (float64); shape: (nBlobsNew,))
                  
                  
        First added:     2013-07-11

        Copyright (C) 2013 Artur Palha
                           pHyFlow
    """
    
    """    
        Reviews:   1- Added option for choosing interpolating kernel. (2013-07-15)
        
                   2- Added xBounds and yBounds in order to have the blobs
                      redistributed to the original grid positions. Otherwise
                      blobs were redistributed to a grid with the same spacing
                      but passing in [0.0, 0.0]. This would introduce a small
                      displacement not desired. (2013-07-18)
    """    

    # determine the number of blobs
    #    nBlobs = xBlob.size    
    
    if xBlob.size < 1: # just in case there are no vortex blobs
        xBlobNew = xBlob
        yBlobNew = yBlob
        wBlobNew = wBlob
    
    else: # if there are vortex blobs
        # compute blob spacing (h)
        # it is assumed that particles are equally spaced in both x and y
        h = sigma * overlap
        
        # each vortex blob is redistributed into a square window of new vortex
        # blobs located at the center of the grid, the size of this window
        # given by the interpolating kernel:
        #   M4' : 4 x 4
        #
        # the following algorithm requires determining the index of the left
        # side of the window for each particle and also the bottom side index
        # of the window for each particle
        #
        # as an example for M4' kernel we have:
        #    
        #       ---------------------------------
        #       |       |       |       |       |
        #       |   X   |   X   |   X   |   X   |
        #       |  (3)  |  (7)  |  (11) |  (15) |
        #       ---------------------------------
        #       |       |       |      O|       |
        #       |   X   |   X   |   X   |   X   |
        #       |  (2)  |  (6)  |  (10) |  (14) |
        #       ---------------------------------
        #       |       |       |       |       |
        #       |   X   |   X   |   X   |   X   |
        #       |  (1)  |  (5)  |  (9)  |  (13) |
        #       ---------------------------------
        #       |       |       |       |       |
        #       |   X   |   X   |   X   |   X   |
        #       |  (0)  |  (4)  |  (8)  |  (12) |
        #       ---------------------------------
        #
        #   0 :: vortex particle location
        #   X :: new vortex particles of the interpolation window
        #
        # the numbering is the internal numbering used when determining the indices
        #
        if interpKernel == blobOptions.M4PRIME_INTERPKERNEL:
            kernelSize = 4
        
        xLIndices = numpy.int64(numpy.ceil((xBlob-xBounds[0]-2.0*h+0.5*h)/h))-1 # the x index
        yLIndices = numpy.int64(numpy.ceil((yBlob-yBounds[0]-2.0*h+0.5*h)/h))-1 # the y index
        
        # convert the indices to numbers greater than 2
        
        #        # start by computing the min and max values of indices
        #        xLIndicesMin = xLIndices.min()
        #        yLIndicesMin = yLIndices.min()
        #        xLIndicesMax = xLIndices.max()
        #        yLIndicesMax = yLIndices.max()

        # generate all particles indices of the window of each particle
        #xKernelIndices = numpy.concatenate((numpy.tile(xLIndices,kernelSize),numpy.tile(xLIndices+1,kernelSize),numpy.tile(xLIndices+2,kernelSize),numpy.tile(xLIndices+3,kernelSize)))
        xKernelIndices = numpy.tile(xLIndices,(kernelSize*kernelSize,1)) + numpy.tile(numpy.arange(0,kernelSize),kernelSize).reshape(kernelSize*kernelSize,1) # slower but general
        #xKernelIndices = xKernelIndices.reshape(kernelSize*kernelSize,nBlobs)        
        yKernelIndices = numpy.tile(yLIndices,(kernelSize*kernelSize,1)) + numpy.repeat(numpy.arange(0,kernelSize),kernelSize).reshape(kernelSize*kernelSize,1) # slower but general

        # clear memory of unecessary data
        del xLIndices
        del yLIndices        
        
        # compute the coordinates of the grid vortex blobs of the interpolating kernel
        #        xKernelCoordinates = xKernelIndices*h;
        #        yKernelCoordinates = yKernelIndices*h;
        
        # compute the weight of the interpolation kernel for each new particle
        kernelValues = _interp_kernel_2d(xKernelIndices + 0.5 + xBounds[0]/h - xBlob/h, yKernelIndices + 0.5 + yBounds[0]/h - yBlob/h, interpKernel, c)
        kernelValues *= numpy.tile(wBlob,(kernelSize*kernelSize,1))
        
        #        # allocate memory space for the sparse matrix
        #        newBlobsMatrix = pysparse.spmatrix.ll_mat(xLIndicesMax - xLIndicesMin + kernelSize + 2,\
        #                                         yLIndicesMax - yLIndicesMin + kernelSize + 2,\
        #                                         kernelSize*kernelSize*nBlobs)
                                         
        # gather data in the sparse matrix
        # repeated vortex blobs are added together
        xKernelIndicesMin = xKernelIndices.min()
        yKernelIndicesMin = yKernelIndices.min()
     
        #        newBlobsMatrix.update_add_at(kernelValues.flatten(),\
        #                           xKernelIndices.flatten() - xKernelIndicesMin + kernelSize/2,\
        #                           yKernelIndices.flatten() - yKernelIndicesMin + kernelSize/2)

        newBlobsMatrix = scipy.sparse.csr_matrix( (kernelValues.flatten(),
                                          (xKernelIndices.flatten() - xKernelIndicesMin + kernelSize/2,\
                                           yKernelIndices.flatten() - yKernelIndicesMin + kernelSize/2))).tocoo()
                           
        # clear memory of unecessary data
        del kernelValues
        del xKernelIndices
        del yKernelIndices
        
        # get the new vortex blob data
        #        wBlobNew,iBlobNew,jBlobNew = newBlobsMatrix.find()
        wBlobNew, iBlobNew, jBlobNew = newBlobsMatrix.data, newBlobsMatrix.row, newBlobsMatrix.col
        
        # clear memory of unecessary data
        del newBlobsMatrix
        
        # convert vortex blob index into coordinates
        xBlobNew = (iBlobNew + xKernelIndicesMin - kernelSize/2)*h + 0.5*h + xBounds[0]
        yBlobNew = (jBlobNew + yKernelIndicesMin - kernelSize/2)*h + 0.5*h + yBounds[0]
                
    
    # return figure handle with the plot                                               
    return xBlobNew,yBlobNew,wBlobNew
    


#def PopulationControl(xBlob,yBlob,wBlob,wThreshold,wRelativeThreshold):
def PopulationControl(xBlob,yBlob,wBlob,wThreshold,wTotalThreshold):
    """
    
        PopulationControl controls the number of particles by discarding
        particles with circulation higher than gThreshold. If total circulation
        of the flagged blobs (flaggedCirculation) is more than 
        wRelativeThreshold then wThreshold is multiplied by 10 and the test is
        performed again.
        
        In the end it is insured that the total circulation did not change more
        than wRelativeThreshold.
        
        Usage
        -----
        
            xBlobNew,yBlobNew,wBlobNew = pop(oldParticles,gThreshold,gRelativeThreshold)
    
        Parameters
        ----------
    
            xBlob :: the x coordinates of the vortex blobs
            -----    (type: numpy.ndarray (float64); shape: (nBlobs,))
                    
            yBlob :: the y coordinates of the vortex blobs
            -----    (type: numpy.ndarray (float64); shape: (nBlobs,))
                     
            wBlob :: the circulations associated to each of the vortex blobs
            -----    (type: numpy.ndarray (float64); shape: (nBlobs,))
                     
    
            wThreshold :: the minimum value of absolute circulation to be
            ----------    considered. Values lower than this are discarded.
                          (type: float64, dimension: single value)
    
            wRelativeThreshold :: the maximum value of the ratio 
            ------------------    flaggedCirculation/totalCirculation that is
                                  acceptable. If
                                  flaggedCirculation/totalCirculation > wRelativeThreshold
                                  wThreshold = wThreshold/10.0 and the process
                                  of flagging blobs is repeated
                                  (type: float64, dimension: single value)

            wTotalThreshold :: the maximum value of the flaggedCirculation 
            ------------------ that is acceptable. If flaggedCirculation  > wRelativeThreshold
                               then wThreshold = wThreshold/10.0 and the process
                               of flagging blobs is repeated
                               (type: float64, dimension: single value)    


        Returns
        -------
            xBlobNew :: the x coordinates of the new vortex blobs
            --------    (type: numpy.ndarray (float64); shape: (nBlobsNew,))
                    
            yBlobNew :: the y coordinates of the new vortex blobs
            --------    (type: numpy.ndarray (float64); shape: (nBlobsNew,))
                     
            wBlobNew :: the circulations associated to each of the new vortex blobs
            --------    (type: numpy.ndarray (float64); shape: (nBlobsNew,))
    
    
        First added:     2013-07-18

        Copyright (C) 2013 Artur Palha
                           pHyFlow
    """
    
    """    
        Reviews:   
    """
    # 1- flag blobs with circulation smaller then wThreshold
    # 2- compute the total circulation of flagged blobs
    # 3- check if total circulation of flagged blobs divided by the total circulation
    #    of all the blobs is smaller than wRelativeThreshold. If so discard
    #    these blobs, if not reduce wThreshold, by dividing it by 10 and go again
    #    to step 1
    
    # compute the total absolute circulation
    #    totalCirculation = numpy.abs(wBlob).sum()
    
    # determine the particles with circulation smaller then gThreshold
    # compute how much absolute circulation is discarded, if 
    # flaggedCirculation/totalCirculation < gRelativeThreshold
    # then discard particles, if not, decrease gThreshold by one order
    # of magnitude
    
    circulationConservationFlag = False # initialize circulation flag that
                                       # determines if the amount of particles
                                       # discarded is acceptable in terms of 
                                       # circulation loss
                                       
    while not circulationConservationFlag:
        flaggedParticles = numpy.abs(wBlob)<wThreshold # determine the blobs to discard
        circulationLoss = (numpy.abs(wBlob[flaggedParticles])).sum() # determine the total circulation loss
        #if (circulationLoss/totalCirculation) < wRelativeThreshold:
        if circulationLoss < wTotalThreshold:
            circulationConservationFlag = True # the circulation change is acceptable
        else:
            wThreshold /=10.0 # decrease the circulation threshold
                                         
    # return the non-discarded blobs
    flaggedParticles = numpy.logical_not(flaggedParticles) # invert the flags, now True means to keep particle
    
    return xBlob[flaggedParticles],yBlob[flaggedParticles],wBlob[flaggedParticles]
    
   
def _M4prime(x,c):
    """
        The one dimensional M4' smooth interpolating kernel given by:
            
                  /
                 | 0                                        if |x| > 2
                 |
        M4'(x) =<  0.5*((2-|x|)^2)*(1-|x|)
                 |      - (c^2)*((2-|x|)^2)*(1-2*|x|)       if |x| \in [1,2]
                 |
                 | 1 - 2.5*(x^2) + 1.5*(|x|^3)
                  \     - (c^2)*(2 - 9*(x^2) + 6*(|x|^3))   if |x| \in [0,1]
                   \    
        
        Usage
        -----
            interpolationWeights = _M4prime(x,c)
            
        Parameters
        ----------
            x :: the points where to compute the interpolation weights
            -    (type: numpy.ndarray (float64); shape: any shape)
            
            c :: the coefficient that models diffusion, given by:
                :math:`c = \\frac{\\sqrt(\\nu * \\Delta T)}{\\Delta x^{2}}`
                 (default value 0.0 (no diffusion))
                 (type: float (64bits); shape: single value)
                    
            
        Returns
        -------
            interpolationWeights :: the interpolation weights
                                   (type: numpy.ndarray (float64); shape: x.shape)
            
                  
        First added:     2013-07-11

        Copyright (C) 2013 Artur Palha
                           pHyFlow
    """
    
    """    
        Reviews:   
    """
    
    # determine the shape of x for reshaping the results in the end
    shape_x = x.shape
    
    # compute the interpolating weights using c++ external function
    # c = 0.0 --> no diffusion
    interpolationWeights = m4kernel_1d(x,c)
    
    # return the interpolation weights reshape to the shape of x because
    # m4kernel_1d return the flattened results of x
    return interpolationWeights.reshape(shape_x)
    
    
    
def _interp_kernel_2d(x,y,interpKernel,c):
    """
        The 2 dimensional interpolating kernels given by tensor product:

            K(x,y) = k(x)*k(y)            
            
        Usage
        -----
            interpolationWeights = _interp_kernel_2d(x,y,interpKernel=0)
            
        Parameters
        ----------
            x :: the x coordinates of the points where to compute the interpolation weights
            -    (type: numpy.ndarray (float64); shape: (n,m))
            
            y :: the y coordinates of the points where to compute the interpolation weights
            -    (type: numpy.ndarray (float64); shape: (n,m))
            
            interpKernel :: the interpolating kernel using for remeshing
            ------------    (type: int (64bits); shape: single value)
            
            c :: the coefficient that models diffusion, given by:
            _    c = \frac{\sqrt(\nu \\times delta)}{\delta x^{2}}
                 (default value 0 (no diffusion))
                 (type: float (64bits); shape: single value)
                    
            
        Returns
        -------
            interpolationWeights :: the interpolation weights
                                   (type: numpy.ndarray (float64); shape: (n,m))
            
                  
        First added:     2013-07-15

        Copyright (C) 2013 Artur Palha
                           pHyFlow
    """
    
    """    
        Reviews:   
    """
    
    # compute the interpolation kernels by tensor product using the chosen
    # interpolation kernel
    
    if interpKernel == blobOptions.M4PRIME_INTERPKERNEL: # M4' interpolating kernel
        interpolationWeights = _M4prime(x,c) * _M4prime(y,c)
    
    # return the interpolation weights
    return interpolationWeights
