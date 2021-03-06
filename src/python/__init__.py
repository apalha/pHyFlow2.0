"""Import vortex and N-S solver modules and options

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
# First added:  2013-05-27                                                                                                          
# Last changed: 2014-03-05
# -*- coding: utf-8 -*-

# load module options
#import options

# load auxiliary modules aux
import aux

# load vortex blob related modules
import blobs

# load panel related modules
import panels

# load lagrangian (blobs + panels) related modules
import lagrangian

# load eulerian (grid solver) modules
import eulerian

# load hybrid related modules
import hybrid

# load input/output modules
import IO


