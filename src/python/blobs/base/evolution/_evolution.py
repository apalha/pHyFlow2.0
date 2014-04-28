"""General routines for computing induced velocities from different types of
vortex kernels.

Several different kernel implementations are implemented in one unified function.
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
# Last changed: 2013-10-09
# -*- coding: utf-8 -*-

__all__ = ['convection','diffusion_wee','diffusion_tutty']

import numpy
import scipy
from pHyFlow.blobs import blobOptions

from pHyFlow.blobs.base.induced import velocity
from pHyFlow.blobs.base.regrid import Regrid


def convection(dt,xBlob,yBlob,wBlob,sigma,k=2,kernel=1,vInf=[0.0,0.0],
               hardware=0,method=1,blocksize=128,integrator=1):
    r"""
        WORK IN PROGRESSS!!!!!!
        
        Compute the new blobs as given by integrating the convection equation:
        
        .. math::
            \frac{dx_{i}}{dt}={\textstyle \sum_{i=1}^{N_{blobs}} v_i\cdot\left( x-x_i \right)}
            
        Usage
        -----
           xBlobNew,yBlobNew,wBlobNew = convection(dt,xBlob,yBlob,wBlob,sigma,k=2,\
                                                   kernel=1,vInf=[0.0,0.0],hardware=0,method=1,\
                                                   blocksize=128,integrator=1)
            
        Parameters
        ----------
            dt :: the size of the time step used in the numerical integration
            --    (type: float (64bits); shape: single value)
            
            xBlob :: the x coordinates of the vortex blobs
            -----    (type: numpy.ndarray (float64); shape: (nBlobs,))
                    
            yBlob :: the y coordinates of the vortex blobs
            -----    (type: numpy.ndarray (float64); shape: (nBlobs,))
                     
            wBlob :: the circulations associated to each of the vortex blobs
            -----    (type: numpy.ndarray (float64); shape: (nBlobs,))
                     
            sigma :: the core size of all the vortex blobs
            -----    (type: float (64bits); shape: single value)
                     
            k :: the core size multiplication constant of all the vortex blobs
                 typical values are 1,2,4 (default value is 2)
                 (type: float (64bits); shape: single value)
                 
            kernel :: the type of kernel of all the vortex blobs
                      available kernels are: 0 ('cutoff'), 1 ('gauss')
                      (default value is gauss kernel)
                      (type: int; shape: single value)
            
            vInf :: the free stream velocity
                    (type: float64; shape: (2,))
                    
            hardware :: the hardware to use to compute the induced velocities
                        can be the 0 (for CPU) or 1 (for GPU)
                        (default value is CPU)
                        (type: int; shape: single value)
            
            method :: the method used to compute the induced velocities can be
                      0 (for FMM) or 1 (for direct calculation)
                      (default value is direct calculation)
                      (type: int; shape: single value)
                                              
            blocksize :: the size of gpu memory block size
                         (default value 128)
                         (type: int; shape: single value)
                         
            integrator :: the numerical integrator to use for integrating the
                          ode
                          available integrator are:
                              0 (forward Euler),
                              1 (Runge-Kutta 4th order)
                     
        Returns
        -------
            xBlobNew :: the x coordinates of the new vortex blobs
            --------    (type: numpy.ndarray (float64); shape: (nBlobs,))
                    
            yBlobNew :: the y coordinates of the new vortex blobs
            --------    (type: numpy.ndarray (float64); shape: (nBlobs,))
            
            wBlobNew :: the circulations associated to each of the new vortex
            --------    blobs
                        (type: numpy.ndarray (float64); shape: (nBlobs,))
            
    :First Added:   2013-05-22
    :Last Modified: 2013-10-10
    :Copyright:     Copyright (C) 2013 Artur Palha, Lento Manickathan **pHyFlow**
    :Licence:       GNU GPL version 3 or any later version 
    
    """
    
    """    
        Reviews:    1) Panel method implemented
    """
    
    if integrator == 0: # forward Euler
        xBlobNew,yBlobNew,wBlobNew = _convection_fe(dt,xBlob,yBlob,wBlob,sigma,
                                                    k,kernel,vInf,
                                                    hardware,method,blocksize)
                                                    
    elif integrator == 1: # Runge-Kutta 4th order
        xBlobNew,yBlobNew,wBlobNew = _convection_rk4(dt,xBlob,yBlob,wBlob,sigma,
                                                     k,kernel,vInf,
                                                     hardware,method,blocksize)
    
    # return the new blobs
    return xBlobNew,yBlobNew,wBlobNew

   
def diffusion_wee(dt,nu,h,xBlob,yBlob,wBlob,sigma,xBounds,yBounds,overlap=0.5,k=2,
                  kernel=1,interpKernel=0,diffusion_method=1):
    """
        Compute the new blobs as given by solving the diffusion equation.
        Diffusion is computed in the following way:
            
               With a special re-distribution process that incorporates
               diffusion, as presented in:
            
                   Wee, D., Ahmed, F., Ghoniem, F., Modified interpolation
                        kernels for treating diffusion and remeshing in vortex
                        methods, Journal of Computational Physics, 213,
                        239--263, 2006.
        
        Usage
        -----
           xBlobNew,yBlobNew,wBlobNew = diffusion_wee(dt,nu,h,xBlob,yBlob,wBlob,
                                                      sigma,
                                                      xBoundsDomain,yBoundsDomain,
                                                      overlap=0.5,k=2,kernel=1)
            
        Parameters
        ----------
            dt :: the size of the time step used in the numerical integration
            --    (type: float (64bits); shape: single value)
            
            nu :: the diffusion coefficient
            --    (type: float (64bits); shape: single value)
            
            h :: the size of the cell associated to the blobs
            -    (type: float (64bits); shape: single value)
            
            xBlob :: the x coordinates of the vortex blobs
            -----    (type: numpy.ndarray (float64); shape: (nBlobs,))
                    
            yBlob :: the y coordinates of the vortex blobs
            -----    (type: numpy.ndarray (float64); shape: (nBlobs,))
                     
            wBlob :: the circulations associated to each of the vortex blobs
            -----    (type: numpy.ndarray (float64); shape: (nBlobs,))
                     
            sigma :: the core size of all the vortex blobs
            -----    (type: float (64bits); shape: single value)
                       
            xBounds :: the x bounds of the domain where the particles were
            -------    distributed originally
                       (type: float (64bits); shape: (2,))
            
            yBounds :: the y bounds of the domain where the particles were
            -------    distributed originally
                       (type: float (64bits); shape: (2,))
            
            overlap :: the overlap between vortex blobs
                       (default value is 0.5) 
                       (type: float (64bits); shape: single value) 
            
            k :: the core size multiplication constant of all the vortex blobs
                 typical values are 1,2,4 (default value is 2)
                 (type: float (64bits); shape: single value)
                 
            kernel :: the type of kernel of all the vortex blobs
                      available kernels are: 0 ('cutoff'), 1 ('gauss')
                      (default value is gauss kernel)
                      (type: int; shape: single value)
                    
        Returns
        -------
            xBlobNew :: the x coordinates of the new vortex blobs
            --------    (type: numpy.ndarray (float64); shape: (nBlobs,))
                    
            yBlobNew :: the y coordinates of the new vortex blobs
            --------    (type: numpy.ndarray (float64); shape: (nBlobs,))
            
            wBlobNew :: the circulations associated to each of the new vortex
            --------    blobs
                        (type: numpy.ndarray (float64); shape: (nBlobs,))
            
    :First Added:   2013-05-22
    :Last Modified: 2014-02-10
    :Copyright:     Copyright (C) 2013 Artur Palha, Lento Manickathan **pHyFlow**
    :Licence:       GNU GPL version 3 or any later version 
    
    """
    
    """    
        Reviews:    1) Change from generic diffusion into diffusion_wee. Each
                       diffusion method will be a specific function.
    """
    
    # Modified interpolation for diffusion
    c = numpy.sqrt(nu*dt)/h

    # Redistribution with Diffusion
    xBlobNew,yBlobNew,wBlobNew = Regrid(xBlob,yBlob,wBlob,sigma,overlap,
                                        xBounds,yBounds,interpKernel,c)
                                       
    # Return the new blobs
    return xBlobNew,yBlobNew,wBlobNew
    

def diffusion_tutty(dt,nu,h,xBlob,yBlob,wBlob,sigma,xBounds,yBounds,overlap=0.5,k=2,
                  kernel=1,interpKernel=0,diffusion_method=1):
    """
        Compute the new blobs as given by solving the diffusion equation.
        Diffusion is computed in the following way:

               With a special re-distribution process that incorporates
               diffusion, as presented in:

                   Tutty, O.R., A simple redistribution vortex method (with accurate body forces),
                   arXiv:1009.0166v1, 2010

        Usage
        -----
           xBlobNew,yBlobNew,wBlobNew = diffusion_tutty(dt,nu,h,xBlob,yBlob,wBlob,
                                                        sigma,
                                                        xBoundsDomain,yBoundsDomain,
                                                        overlap=0.5,k=2,kernel=1)

        Parameters
        ----------
            dt :: the size of the time step used in the numerical integration
            --    (type: float (64bits); shape: single value)

            nu :: the diffusion coefficient
            --    (type: float (64bits); shape: single value)

            h :: the size of the cell associated to the blobs
            -    (type: float (64bits); shape: single value)

            xBlob :: the x coordinates of the vortex blobs
            -----    (type: numpy.ndarray (float64); shape: (nBlobs,))

            yBlob :: the y coordinates of the vortex blobs
            -----    (type: numpy.ndarray (float64); shape: (nBlobs,))

            wBlob :: the circulations associated to each of the vortex blobs
            -----    (type: numpy.ndarray (float64); shape: (nBlobs,))

            sigma :: the core size of all the vortex blobs
            -----    (type: float (64bits); shape: single value)

            xBounds :: the x bounds of the domain where the particles were
            -------    distributed originally
                       (type: float (64bits); shape: (2,))

            yBounds :: the y bounds of the domain where the particles were
            -------    distributed originally
                       (type: float (64bits); shape: (2,))

            overlap :: the overlap between vortex blobs
                       (default value is 0.5)
                       (type: float (64bits); shape: single value)

            k :: the core size multiplication constant of all the vortex blobs
                 typical values are 1,2,4 (default value is 2)
                 (type: float (64bits); shape: single value)

            kernel :: the type of kernel of all the vortex blobs
                      available kernels are: 0 ('cutoff'), 1 ('gauss')
                      (default value is gauss kernel)
                      (type: int; shape: single value)

        Returns
        -------
            xBlobNew :: the x coordinates of the new vortex blobs
            --------    (type: numpy.ndarray (float64); shape: (nBlobs,))

            yBlobNew :: the y coordinates of the new vortex blobs
            --------    (type: numpy.ndarray (float64); shape: (nBlobs,))

            wBlobNew :: the circulations associated to each of the new vortex
            --------    blobs
                        (type: numpy.ndarray (float64); shape: (nBlobs,))

    :First Added:   2014-04-14
    :Last Modified: 2014-04-14
    :Copyright:     Copyright (C) 2014 Artur Palha, Lento Manickathan **pHyFlow**
    :Licence:       GNU GPL version 3 or any later version

    """

    """
        Reviews:

    """

    if xBlob.size < 1: # just in case there are no vortex blobs
        xBlobNew = xBlob
        yBlobNew = yBlob
        wBlobNew = wBlob

    else: # if there are vortex blobs
        # compute blob spacing (h)
        # it is assumed that particles are equally spaced in both x and y
        h = sigma * overlap


        # Generate the coordinates of the kernel nodes

        # Determine the x,y lower index
        xLIndices = numpy.int64(numpy.ceil((xBlob-xBounds[0]-2.0*h+0.5*h)/h))-1 # the x index
        yLIndices = numpy.int64(numpy.ceil((yBlob-yBounds[0]-2.0*h+0.5*h)/h))-1 # the y index

        # Determine the reference index : xi
        xIndices = xLIndices + 1 # the x index
        yIndices = yLIndices + 1 # the y index

        # convert the indices to numbers greater than 2

        # generate all particles indices of the window of each particle
        #xKernelIndices = numpy.concatenate((numpy.tile(xLIndices,kernelSize),numpy.tile(xLIndices+1,kernelSize),numpy.tile(xLIndices+2,kernelSize),numpy.tile(xLIndices+3,kernelSize)))
        xKernelIndices = numpy.tile(xLIndices,(4*4,1)) + numpy.tile(numpy.arange(0,4),4).reshape(4*4,1) # slower but general
        #xKernelIndices = xKernelIndices.reshape(kernelSize*kernelSize,nBlobs)
        yKernelIndices = numpy.tile(yLIndices,(4*4,1)) + numpy.repeat(numpy.arange(0,4),4).reshape(4*4,1) # slower but general

        # clear memory of unnecessary data
        del xLIndices
        del yLIndices


        # Description of the Tutty redistribution method (with diffusion)
        # --------------------------------------------------------------
        #
        # - particles are placed at x = xNu.
        # - Delta = (xNu - xi)/h
        # - xi = i*h (in our case == > xi = i*h + 0.5*h)


        # compute the characteristic diffusion distance
        hNu = numpy.sqrt(dt*nu)

        # Determine the coordinates of the particle reference point i
        xICoordinates = xIndices*h + 0.5*h + xBounds[0]
        yICoordinates = yIndices*h + 0.5*h + yBounds[0]

        del xIndices
        del yIndices

        # Determine the Delta of each particles  # !!!!!! MAKE SURE THIS IS CORRECT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        # Note: xNu = xBlob, yNu = yBlob
        xDelta = (xBlob - xICoordinates)/h
        yDelta = (yBlob - yICoordinates)/h

        # Determine the Delta1 of each particles
        xDelta1 = xDelta - 1.0
        yDelta1 = yDelta - 1.0

        # Determine the f_i and g_{i+1} redistribution equations - x-direction
        fx_i   = 1.0 - 2.0*(hNu/h)**2.0 - xDelta**2.0
        gx_ip1 = 1.0 - 2.0*(hNu/h)**2.0 - xDelta1**2.0

        # Determine the f_j and g_{j+1} redistribution equations - y-direction
        fy_i   = 1.0 - 2.0*(hNu/h)**2.0 - yDelta**2.0
        gy_ip1 = 1.0 - 2.0*(hNu/h)**2.0 - yDelta1**2.0


        # Perform linear combinations of the two basis solutions. - x-direction
        # F_k = (1 - Delta)*f_k + Delta*g_k,  where k = i-1,...,i+2
        F_im1 = (1.0-xDelta) * (0.5 * (1.0-fx_i-xDelta))                                          # F_{k=i-1}
        F_i   = (1.0-xDelta) * fx_i                     + xDelta * (0.5 * (1.0-gx_ip1-xDelta1))   # F_{k=i}
        F_ip1 = (1.0-xDelta) * (0.5 * (1.0-fx_i+xDelta))  + xDelta * gx_ip1                       # F_{k=i+1}
        F_ip2 =                                        xDelta * (0.5 * (1.0-gx_ip1+xDelta1))      # F_{k=i+2}

        # Perform linear combinations of the two basis solutions. - x-direction
        # G_k = (1 - Delta)*f_k + Delta*g_k,  where k = j-1,...,j+2
        G_im1 = (1.0-yDelta) * (0.5 * (1.0-fy_i-yDelta))                                          # F_{k=j-1}
        G_i   = (1.0-yDelta) * fy_i                     + yDelta * (0.5 * (1.0-gy_ip1-yDelta1))   # F_{k=j}
        G_ip1 = (1.0-yDelta) * (0.5 * (1.0-fy_i+yDelta))  + yDelta * gy_ip1                       # F_{k=j+1}
        G_ip2 =                                         yDelta * (0.5 * (1.0-gy_ip1+yDelta1))     # F_{k=j+2}

        # Stack and tile the 1-D redistribution weights for the 4x4 2-D kernel - x-direction
        FWeights = numpy.tile(numpy.vstack((F_im1,F_i,F_ip1,F_ip2)),(4,1)).flatten()
        # Stack and tile the 1-D redistribution weights for the 4x4 2-D kernel - y-direction
        GWeights = numpy.vstack((numpy.tile(G_im1,(4,1)),
                                 numpy.tile(G_i,  (4,1)),
                                 numpy.tile(G_ip1,(4,1)),
                                 numpy.tile(G_ip2, (4,1)))).flatten()






        # compute the coordinates of the grid vortex blobs of the interpolating kernel
        #        xKernelCoordinates = xKernelIndices*h;
        #        yKernelCoordinates = yKernelIndices*h;

        # compute the weight of the interpolation kernel for each new particle
        kernelValues = FWeights*GWeights
        kernelValues *= numpy.tile(wBlob,(4*4,1)).flatten()

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
                                          (xKernelIndices.flatten() - xKernelIndicesMin + 4/2,\
                                           yKernelIndices.flatten() - yKernelIndicesMin + 4/2))).tocoo()

        # clear memory of unnecessary data
        del kernelValues
        del xKernelIndices
        del yKernelIndices

        # get the new vortex blob data
        #        wBlobNew,iBlobNew,jBlobNew = newBlobsMatrix.find()
        wBlobNew, iBlobNew, jBlobNew = newBlobsMatrix.data, newBlobsMatrix.row, newBlobsMatrix.col

        # clear memory of unecessary data
        del newBlobsMatrix

        # convert vortex blob index into coordinates
        xBlobNew = (iBlobNew + xKernelIndicesMin - 4/2)*h + 0.5*h + xBounds[0]
        yBlobNew = (jBlobNew + yKernelIndicesMin - 4/2)*h + 0.5*h + yBounds[0]


    # return figure handle with the plot
    return xBlobNew,yBlobNew,wBlobNew


def _convection_fe(dt,xBlob,yBlob,wBlob,sigma,k,kernel,vInf,hardware,method,
                   blocksize):
                   
    r"""
        Compute the new blobs as given by integrating the convection equation:
        
        .. math::
            \frac{dx_{i}}{dt}={\textstyle \sum_{i=1}^{N_{blobs}} v_i\cdot\left( x-x_i \right)}
        
        using the forward Euler method.
        
        Usage
        -----
           xBlobNew,yBlobNew,wBlobNew = convection_fe(dt,xBlob,yBlob,wBlob,sigma,k,
                                                   kernel,vInf,hardware,method,blocksize)
            
        Parameters
        ----------
            dt :: the size of the time step used in the numerical integration
            --    (type: float (64bits); shape: single value)
            
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
                      
            vInf :: the free stream velocity
            ----    (type: float64; shape: (2,))
            
            hardware :: the hardware to use to compute the induced velocities
            --------    can be the 0 (for CPU) or 1 (for GPU)
                        (default value is CPU)
                        (type: int; shape: single value)
            
            method :: the method used to compute the induced velocities can be
            ------    0 (for FMM) or 1 (for direct calculation)
                      (default value is direct calculation)
                      (type: int; shape: single value)
                                              
            blocksize :: the size of gpu memory block size
            ---------    (default value 128)
                         (type: int; shape: single value)
                         
            integrator :: the numerical integrator to use for integrating the
            ---------     ode
                          available integrator are: 0 (forward Euler),
                          1 (Runge-Kutta 4th order)
                     
        Returns
        -------
            xBlobNew :: the x coordinates of the new vortex blobs
            --------    (type: numpy.ndarray (float64); shape: (nBlobs,))
                    
            yBlobNew :: the y coordinates of the new vortex blobs
            --------    (type: numpy.ndarray (float64); shape: (nBlobs,))
            
            wBlobNew :: the circulations associated to each of the new vortex
            --------    blobs
                        (type: numpy.ndarray (float64); shape: (nBlobs,))
            
    :First Added:   2013-05-22
    :Last Modified: 2013-10-10
    :Copyright:     Copyright (C) 2013 Artur Palha, Lento Manickathan **pHyFlow**
    :Licence:       GNU GPL version 3 or any later version 

    """
    
    """    
        Reviews:       1) Panel method implemented
    """
          
    # Vortex method WITHOUT solid boundaries
    print "convection : vortex methods WITHOUT panels."   
   
    # compute induced velocities
    # use the options the user gave as input and compute the velocities at
    # blob centers
    vx,vy = velocity(xBlob,yBlob,wBlob,sigma,k=k,kernel=kernel,\
             xEval=None,yEval=None,hardware=hardware,blocksize=blocksize,method=method)
      
    # advance blobs in time with forward Euler
    # x_new = x_old + dt * vx
    xBlobNew = xBlob + dt * (vx+vInf[0]);
    yBlobNew = yBlob + dt * (vy+vInf[1]);
                                            

    # circulations are conserved in the convection step
    wBlobNew = wBlob
    
    return xBlobNew,yBlobNew,wBlobNew
             


def _convection_rk4(dt,xBlob,yBlob,wBlob,sigma,k,kernel,
                    vInf,hardware,method,blocksize):
    r"""
        Compute the new blobs as given by integrating the convection equation:
        
        .. math::
            \frac{dx_{i}}{dt}={\textstyle \sum_{i=1}^{N_{blobs}} v_i\cdot\left( x-x_i \right)}
            
        using the Runge-Kutta of 4th order method.
        
        Usage
        -----
           xBlobNew,yBlobNew,wBlobNew = convection_rk4(dt,xBlob,yBlob,wBlob,sigma,k,
                                                   kernel,vInf,hardware,method,blocksize)
            
        Parameters
        ----------
            dt :: the size of the time step used in the numerical integration
            --    (type: float (64bits); shape: single value)
            
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
                      
            vInf :: the free stream velocity
            ----    (type: float64; shape: (2,))
            
            hardware :: the hardware to use to compute the induced velocities
            --------    can be the 0 (for CPU) or 1 (for GPU)
                        (default value is CPU)
                        (type: int; shape: single value)
            
            method :: the method used to compute the induced velocities can be
            ------    0 (for FMM) or 1 (for direct calculation)
                      (default value is direct calculation)
                      (type: int; shape: single value)
                                              
            blocksize :: the size of gpu memory block size
            ---------    (default value 128)
                         (type: int; shape: single value)
                         
            integrator :: the numerical integrator to use for integrating the
            ---------     ode
                          available integrator are: 0 (forward Euler),
                          1 (Runge-Kutta 4th order)
                     
        Returns
        -------
            xBlobNew :: the x coordinates of the new vortex blobs
            --------    (type: numpy.ndarray (float64); shape: (nBlobs,))
                    
            yBlobNew :: the y coordinates of the new vortex blobs
            --------    (type: numpy.ndarray (float64); shape: (nBlobs,))
            
            wBlobNew :: the circulations associated to each of the new vortex
            --------    blobs
                        (type: numpy.ndarray (float64); shape: (nBlobs,))
            
        First added:     2013-09-29

        Copyright (C) 2013 Artur Palha
                           pHyFlow
    """
    
    """    
        Reviews:    1) panel method implement to enforce no-slip b.c at wall
    """
    
            
    # advance blobs in time with RK4
    # x_new = x_old + (1/6) * dt * (k1 + 2*k2 + 2*k3 + k4)   (1)
    #
    # k1 = v(t,x_old)                                        (2)
    # k2 = v(t+0.5*dt,x_old + 0.5*dt*k1)                     (3)
    # k3 = v(t+0.5*dt,x_old + 0.5*dt*k2)                     (4)
    # k4 = v(t+dt,x_old + dt*k3)                             (5)

    # in order to compute more efficiently the new positions, expression (1)
    # will be computed in steps:
    #
    # x_new_1 = k1                          (step 1)
    # x_new_2 = x_new_1 + 2*k2              (step 2)
    # x_new_3 = x_new_2 + 2*k3              (step 3)
    # x_new_4 = x_new_3 + k4                (step 4)
    # x_new_5 = (1/6)*dt*x_new_4 + x_old    (step 5)
        
    # allocate memory space for new blobs
    xBlobNew = numpy.zeros(xBlob.shape)
    yBlobNew = numpy.zeros(yBlob.shape)
    
    # STEP 1 :: x_new_1 = k1
    
    # compute kx1 and ky1 using (2)
    vx,vy = velocity(xBlob,yBlob,wBlob,sigma,k=k,kernel=kernel,\
             xEval=None,yEval=None,hardware=hardware,blocksize=blocksize,method=method) # evaluate the velocities at (xBlob,yBlob)
    kxTemp = vx + vInf[0] # compute kx1 using (2)
    kyTemp = vy + vInf[1] # compute ky1 using (2)

    # update the positions using (step 1)
    xBlobNew += kxTemp # kxTemp is now kx1
    yBlobNew += kyTemp # kyTemp is now ky1
    

    # STEP 2 :: x_new_2 = x_new_1 + 2*k2
    
    # compute kx2 and kx2 using (3)
    vx,vy = velocity(xBlob + 0.5*dt*kxTemp,yBlob + 0.5*dt*kyTemp,wBlob,sigma,k=k,kernel=kernel,\
             xEval=None,yEval=None,hardware=hardware,blocksize=blocksize,method=method) # evaluate the velocities at (xBlob+0.5*dt*kx1,yBlob+0.5*dt*ky1)
                                                                                        # note that kx1 = kxTemp and ky1 = kyTemp at this stage
    kxTemp = vx + vInf[0] # compute kx2 using (3)
    kyTemp = vy + vInf[1] # compute ky2 using (3)
    
    # uptade the positions using (step 2)
    xBlobNew += 2.0*kxTemp # kxTemp is now kx2
    yBlobNew += 2.0*kyTemp # kyTemp is now ky2


    # STEP 3 :: x_new_3 = x_new_2 + 2*k3
    
    # compute kx3 and kx3 using (4)
    vx,vy = velocity(xBlob + 0.5*dt*kxTemp,yBlob + 0.5*dt*kyTemp,wBlob,sigma,k=k,kernel=kernel,\
             xEval=None,yEval=None,hardware=hardware,blocksize=blocksize,method=method) # evaluate the velocities at (xBlob+0.5*dt*kx1,yBlob+0.5*dt*ky1)
                                                                                        # note that kx2 = kxTemp and ky2 = kyTemp at this stage
    kxTemp = vx + vInf[0] # compute kx3 using (4)
    kyTemp = vy + vInf[1] # compute ky3 using (4)
    
    # uptade the positions using (step 3)
    xBlobNew += 2.0*kxTemp # kxTemp is now kx3
    yBlobNew += 2.0*kyTemp # kyTemp is now ky3


    # STEP 4 :: x_new_4 = x_new_3 + k4
    
    # compute kx4 and kx4 using (5)
    vx,vy = velocity(xBlob + dt*kxTemp,yBlob + dt*kyTemp,wBlob,sigma,k=k,kernel=kernel,\
             xEval=None,yEval=None,hardware=hardware,blocksize=blocksize,method=method) # evaluate the velocities at (xBlob+0.5*dt*kx1,yBlob+0.5*dt*ky1)
                                                                                        # note that kx3 = kxTemp and ky3 = kyTemp at this stage
    kxTemp = vx + vInf[0] # compute kx3 using (5)
    kyTemp = vy + vInf[1] # compute ky3 using (5)
    
    # uptade the positions using (step 3)
    xBlobNew += kxTemp # kxTemp is now kx4
    yBlobNew += kyTemp # kyTemp is now ky4


    # STEP 5 :: x_new_5 = (1/6)*dt*x_new_4 + x_old
    
    # x_new_5a = (1/6)*dt*x_new_4
    xBlobNew *= (1.0/6.0) * dt
    yBlobNew *= (1.0/6.0) * dt
    
    # x_new_5b = x_new_5a + x_old
    xBlobNew += xBlob 
    yBlobNew += yBlob
        
        
    # circulations are conserved in the convection step
    wBlobNew = wBlob.copy()
    
    # Return new blob positions
    return xBlobNew,yBlobNew,wBlobNew
             
