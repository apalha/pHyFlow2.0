"""This module imports all Gauss kernels.
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

# induced velocity
from _vv2parGauss_C_cpu import vv2parGauss_C_cpu
from vv2parGauss_Cython_cpu import vv2parGauss_Cython_cpu
from _vv2parGauss_gpu import vv2parGauss_gpu
from _fmm2d_py_cpu import fmm2d_py as fmm2d_cpu
from _fmm2d_py_gpu import fmm2d_py as fmm2d_gpu

# induced vorticity
from _ww2parGauss_C_cpu import ww2parGauss_C_cpu
from _ww2parGauss_gpu import ww2parGauss_gpu
