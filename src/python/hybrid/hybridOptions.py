#-*- coding: utf-8 -*-
r"""
Define hybrid options
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
# First added:  2014-02-21                                                                                                          
# Last changed: 2014-02-28
# -*- coding: utf-8 -*-
##############################################################################

"""    

"""

# When should we correct the lagrangian domain.
# Stock _[1] and Daeninck _[2] correct the lagrangian field at the start
#
# .. [1] Stock, M., Gharakhani, A., & Stone, C. (2010). Modeling Rotor Wakes 
#        with a Hybrid OVERFLOW-Vortex Method on a GPU Cluster. AIAA Applied 
#        Aerodynamics
# .. [2] Daeninck, G. (2006). D EVELOPMENTS IN HYBRID APPROACHES Vortex method 
#        with known separation location Vortex method with near-wall Eulerian 
#        solver RANS-LES coupling. Universit´ catholique de Louvain.
#

ADJUST_LAGRANGIAN_AT = {'default': 'start',
                       'available': ('start', 'end')}
                       
# We could also not couple the problem, but just evolve both solutions
ADJUST_LAGRANGIAN = {'default': True,
                     'available': (True, False)}
                    
                   
# Determine the eulerian initian conditions
EULERIAN_INITIAL_CONDITIONS = {'default': 'lagrangian_field',
                              'available': ('lagrangian_field','eulerian_field') }


# Interpolation parameters
INTERPOLATION_ALGORITHM = {'default': 'scipy_griddata',
                           'available': ('scipy_griddata','structured_probes')}
                      
# Method of interpolation                      
INTERPOLATION_METHOD = {'default': 'linear',
                        'available': ('linear')}                   