r"""
Define vortex-panel options
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
# Last changed: 2014-02-21
# -*- coding: utf-8 -*-
##############################################################################

"""    

"""

# Panel strength update statement
# 'constant': the panel strength is only calculated at the beginning of the 
#             convection step. It is kept constant throughout the RK4 time
#             integration
# 'varying' : the panel strength is recalculated during each sub-step of 
#             rk4 time integration
PANEL_STRENGTH_UPDATE = {'default': 'constant',
                         'available': ('constant', 'varying')}
