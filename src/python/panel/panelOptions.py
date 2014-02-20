r"""
Define panel options
"""
##############################################################################
# Copyright (C) 2014 Lento Manickathan
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
# First added:  2014-02-20                                                                                                          
# Last changed: 2014-02-17
# -*- coding: utf-8 -*-
##############################################################################

"""    

"""

import pHyFlow.options as _options

# Panel kernels
PANEL_KERNEL = {'default'   : 'csv', # Constant Strength Vortex panel
                'available' : ('csv', 'css')} # Constant Strength Vortex, Constant Strength Source Panels

# Problem types
PROBLEM_TYPE = {'default'   : 'moving',
                'available' : ('fixed', 'moving')}

# Velocity computation parameters
VELOCITY_COMPUTATION_METHOD = {'default'    : 'direct',
                               'available'  : ('direct')}

# Velocity computation hardwares
VELOCITY_COMPUTATION_HARDWARE = {'default'      : 'cpu',
                                 'available'    : ('cpu')}

# Solver computation parameters
SOLVER_COMPUTATION_METHOD = {'default'  : 'bicgstab',
                             'available': ('direct', 'gmres', 'bicgstab')}

# Solver computation tolerance
SOLVER_COMPUTATION_TOL = {'default': 1e-12}

SOLVER_COMPUTATION_ASSEMBLE = {'default'  : 'all',
                               'available': ('all')}


#---------------------------------------------------------------------------
# Boolean options list

# Panel kernels
_PANEL_KERNEL_CSV = 1
_PANEL_KERNEL_CSS = 0

# Problem types
_PROBLEM_TYPE_FIXED = 0
_PROBLEM_TYPE_MOVING = 1

# Velocity computation
_VELOCITY_COMPUTATION_METHOD_DIRECT = _options.DIRECT_METHOD
_VELOCITY_COMPUTATION_HARDWARE_CPU  = _options.CPU_HARDWARE





