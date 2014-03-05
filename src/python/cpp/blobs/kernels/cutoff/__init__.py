"""This module imports all cutoff kernels.
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
# Last changed: 2013-05-22

# induced velocities
from _vv2par_C_cpu import vv2par_C_cpu
from vv2par_Cython_cpu import vv2par_Cython_cpu
from _vv2par_gpu import vv2par_gpu

# induced vorticites
from _ww2par_C_cpu import ww2par_C_cpu
from _ww2par_gpu import ww2par_gpu
