"""This module installs Python source files, compiling all necessary builtin functions
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

import os
from distutils import dir_util


#############################################################################
# Install to libPath

# compile and copy optimized external functions
os.chdir(os.path.join(os.curdir,externCodeDirName))
execfile('Make.py')

# copy Python source files and compiled externCode to libRoot
os.chdir(srcSourcePath)
dir_util.copy_tree('python',libRoot)

# copy example files to libExamples
os.chdir(srcSourcePath)
dir_util.copy_tree('examples',libExamples)

#############################################################################
