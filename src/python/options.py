"""Define package wide options

Several different kernel implementations are implemented in one unified function.
"""
##############################################################################
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
# First added:  2013-05-27                                                                                                          
# Last changed: 2013-06-20
# -*- coding: utf-8 -*-
##############################################################################

"""    
    Reviews:    1. Added FMM_METHOD and DIRECT_METHOD options 
                   for computing the induced velocities. Added
                   FE_INTEGRATOR and RK4_INTEGRATOR options for time integration.
                   (2013-05-28)
                
                2. Added W_PLOT (vorticity) and V_PLOT (velocity)
                   options for plotting the induced fields. Added
                   WFUNCTION_PLOT and WBLOB_PLOT options for plotting vorticity
                   field. (2013-06-27)
"""

# general variables to simplify choosing options

# vortex kernel options
GAUSS_KERNEL = 1 # the gaussian kernel
CUTOFF_KERNEL = 0 # the cutoff kernel

# hardware options
CPU_HARDWARE = 0 # the cpu hardware option
GPU_HARDWARE = 1 # the gpu hardware option

# method options for induced velocity
FMM_METHOD = 0 # the FMM method (N.log(N) complexity) option
DIRECT_METHOD = 1 # the direct calculation (N.N complexity) option

# integrator options for computing the evolution of the blobs
FE_INTEGRATOR = 0 # forward Euler (for testing only, highly not adviseable)
RK4_INTEGRATOR = 1 # Runge-Kutta 4th order

# plot options to choose induced field
W_FIELD = 0 # vorticity field
V_FIELD = 1 # velocity field

# plot options for vorticity field
WBLOB_PLOT = 0 # use blob circulation to compute vorticity field
WFUNCTION_PLOT = 1 # use blob functions to compute the vorticity field

# compute options for vorticity field
WBLOB_VORTICITY = 0 # use blob circulation to compute vorticity field
WFUNCTION_VORTICITY = 1 # use blob functions to compute the vorticity field

# define the type of plot
PYLAB_PLOT = 0 # pylab plot
DOLFIN_PLOT = 1 # dolfin plot

# define the type of interpolation kernel
M4PRIME_INTERPKERNEL = 0

# define the type of diffusion method
REGRID_WEE_DIFFUSION = 0

# package wide runtime options
vortex = {}

# hardware
#    cpu :: all routines use cpu version if available
#    gpu :: all routines use gpu version if available
vortex['hardware'] = {'value':0,'description':str('-----------------------------------------------------------------------\n  \'hardware\'      | %d :: all routines use cpu version if available\n                  |\n                  | %d :: all routines use gpu version if available' % (CPU_HARDWARE, GPU_HARDWARE))} 
# gpublocksize
#    the size of the gpu block
vortex['gpublocksize'] = {'value':256,'description':'-----------------------------------------------------------------------\n  \'gpublocksize\'  | [int] :: an integer specifying the size of the\n                  |          memory block size on the gpu             '} 


# verbose options

# time integrators
time_integrator_options = ('euler','rk4')

# hardware
hardware_options = ('cpu','gpu')

# biot savart methods
biot_savart_options = ('direct','fmm')

# diffusion method options
diffusion_method_options = ('regrid_wee')
                           

def help_vortex():
    """
        Print the explanation of all options.
        
        Usage
        -----
            help()
            
        Parameters
        ----------
            None

        Returns
        -------
            String with the description of all options
            
        Copyright (C) 2013 Artur Palha
                           pHyFlow
    """
    
    print '-----------------------------------------------------------------------'
    print ' Parameter Name   | Options                                        '
    for optionkey in vortex.keys():
        print vortex[optionkey]['description']
        
    print '-----------------------------------------------------------------------' 
                               
                           