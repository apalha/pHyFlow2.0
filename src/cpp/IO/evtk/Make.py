"""Python Makefile to automate the build and installation of evtk
"""
# Copyright (C) 2014 Artur Palha                                                                                                     
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
# First added:  2014-02-19                                                                                                          
# Last changed: 2013-02-19
# -*- coding: utf-8 -*-

import os
import shutil

# build evtk
os.system('python setup.py build')

# copy evtk build directory into the python/IO directory

# first find the name of the directory with the compiled evtk
directoriesBuild = os.listdir(os.path.join(os.curdir,'build'))
# loop over the directories and find the one starting with lib
for evtkDir in directoriesBuild:
    if evtkDir[0:3] == 'lib':
        break

evtkBuildPath = os.path.join(os.curdir,'build',evtkDir,'evtk')

# then copy
if os.path.isdir(os.path.join(libIOPath,'__evtk')): # if the directory __evtk already exists delete it
    shutil.rmtree(os.path.join(libIOPath,'__evtk'))

shutil.copytree(evtkBuildPath,os.path.join(libIOPath,'__evtk'))


