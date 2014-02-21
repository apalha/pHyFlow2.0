"""This module installs pHyFlow, compiling all necessary builtin functions
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

import commands
import os
import sys
import shutil
from distutils import dir_util

#############################################################################
# User defined installation parameters

moduleName = 'pHyFlow' # name of the module, used for root folder

libParent = '/media/DATAPART1/programs/lib/python' # root install directory

#############################################################################


#############################################################################
# Determine the inputs to the install process

# check the inputs
input_read = True
input_n = 1

while input_read and len(sys.argv)>1: # if no inputs are given, skip input read
    if sys.argv[input_n] == 'clean': # delete any previously installed version (no installation is performed)
        cleanFlag = True
        buildFlag = False
        input_n += 1    # increment input index
        
    elif sys.argv[input_n] == 'build': # delete any previously installed versions and rebuild all external functions
        cleanFlag = True        
        buildFlag = True
        input_n += 1    # increment input index
        
    elif sys.argv[input_n] == 'install': # install pHyFlow
        cleanFlag = False        
        buildFlag = True
        input_n += 1    # increment input index
        
    elif sys.argv[input_n] == 'help': # print the help
        print 'usage: python install.py [option1] ... [optionN]'
        print '       by default if any previous version is installed newer versions are overwritten but not rebuilt, except if source files have changed.'
        print '\nOptions:'
        print '   clean     : delete any previously installated version (no installation is performed) '
        print '   build     : delete any previously installated version, rebuild all external functions and install pHyFlow'
        print '   install   : install pHyFlow, external functions are only rebuilt if they haven\'t been built yet'
        input_n += 1    # increment input index
        
    else:
        print str(sys.argv[input_n])
        input_n += 1 # increment two indices since this input has two items
        
    if input_n >= len(sys.argv): # if all the inputs have been scanned
        input_read = False   # stop reading them

if len(sys.argv) == 1: # if no input is given just do install
    cleanFlag = False
    buildFlag = True
    
#############################################################################


#############################################################################
# Automated installation parameters

# root path (where install process starts)
startInstallPath = os.path.abspath(os.curdir) # get the absolute path of root install directory

# install paths
libRoot = os.path.join(libParent,moduleName) # define the library root directory

rootSourcePath = os.path.join(libRoot,'build') # get the absolute path of root build directory

srcSourcePath = os.path.join(rootSourcePath,'src') # get the absolute path of the source install path 

pythonSrcSourcePath = os.path.join(srcSourcePath,'python') # get the absolute path of the python install path 
                                                           # where all compiled vortex functions are stored
libExamples = os.path.join(libRoot,'examples') # define the examples root directory                                

libVortexPath = os.path.join(libRoot,'vortex') # define the vortex library root directory

libIOPath = os.path.join(libRoot,'IO') # define the IO library root directory

externCodeDirName = 'cpp' # define the directory where to store all compiled external libraries

externCodeVortexPath = os.path.join(pythonSrcSourcePath,externCodeDirName,'vortex') # define the vortex library compiled external libraries root directory
                                                             
externCodePanelPath = os.path.join(pythonSrcSourcePath,externCodeDirName,'panel') # define the panel library compiled eternal 

#############################################################################


#############################################################################
# Set up directory structure

print '--- pHyFlow installation ---------------------------------------------'
print ' '
print ' install path :: ' + libRoot
print ' '
print '----------------------------------------------------------------------'

# ask the user to proceed with installation or not
gatekeeper = True
while gatekeeper:
    key_pressed = raw_input('proceed to install? y/n > ')
    if key_pressed == 'y' or key_pressed == 'Y':
        gatekeeper = False
    elif key_pressed == 'n' or key_pressed == 'N':
        sys.exit('User terminated installation')
    else:
        print '\n you must press either y or n!\n'
        
# check if parent directory of root directory already exists if not exit
if not os.path.isdir(libParent):
    raise BaseException('\n ' + moduleName + ' library root dir does not exist :: ' + libParent + '\n create this directory and try again')
    
# check if root directory already exists if not create it
if not os.path.isdir(libRoot):
    os.mkdir(libRoot)

# check if examples directory already exists if not create it
if not os.path.isdir(libExamples):
    os.mkdir(libExamples)    

# check if build directory already exists if not create it
if not os.path.isdir(rootSourcePath):
    os.mkdir(rootSourcePath)    

print '\n Installation process started\n'

# copy all the files to the build directory
dir_util.copy_tree(startInstallPath,rootSourcePath)

# execute src install file
os.chdir(os.path.join(rootSourcePath,'src'))
execfile('install.py')

#############################################################################
