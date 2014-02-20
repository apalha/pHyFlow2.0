r"""
Define navier-stokes options
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
# Last changed: 2014-02-20
# -*- coding: utf-8 -*-
##############################################################################

"""    

"""

# Fluid/Boundary Domain integer ID
ID_FLUID_DOMAIN         = 1 # the fluid domain
ID_NOSLIP_BOUNDARY      = 2 # no-slip boundary
ID_EXTERNAL_BOUNDARY    = 3 # external boundary (where dirichlet b.c is applied)
ID_PRESSURE_OUTLET      = 4 # pressure outlet

# NS Solvers
SOLVER = {'default':'ipcs',             # Default solver
          'solvers': ('chorin','ipcs')} # Available solvers

# Partial differential equation solver parameters                            
FORM_COMPILER = {'cpp_optimize': True}

KRYLOV_SOLVER = {'absolute_tolerance':  1e-25,
                 'relative_tolerance':  1e-12,
                 'monitor_convergence': False}

ALLOW_EXTRAPOLATION = True

SET_LOG_ACTIVE = False 

