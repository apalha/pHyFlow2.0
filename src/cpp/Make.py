"""Python Makefile to automate the build and installation of all C Python modules
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
# -*- coding: utf-8 -*-

import os

# compile and install all external functions

# Start by getting of list of the contents of this directory
dirNames = os.listdir(os.curdir)

# Loop over all the contents
for dirName in dirNames:
    # But restrict to the ones that are directories
    if os.path.isdir(dirName):
        # change to that directory and execute the Make.py file inside it
        os.chdir(os.path.join(os.curdir,dirName))
        execfile('Make.py')
        # change back to the cpp directory
        os.chdir(os.path.join(rootSourcePath,'src',externCodeDirName))