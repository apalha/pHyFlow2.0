"""Python Makefile to automate the build and installation of vortex C Python modules
"""
# Copyright (C) 2013 Artur Palha, Lento Manickathan                                                                                                    
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
# Last changed: 2013-02-12
# -*- coding: utf-8 -*-

import commands
import os
import sys
import shutil

srcPath = os.path.abspath(os.path.curdir)

kernelsPath = os.path.join(externCodePanelsPath, 'kernels')


# check if libPath already exists if not exit
if not os.path.isdir(libRoot):
    raise BaseException('library path does not exist')

# update the panel method functions (Cython)
panelKernels = [{'name'         :'Vortex Panel Kernel', 
                 'source_type'  :'Cython',
                 'source_path'  : os.path.join(srcPath,'cython','vortexPanel_Cython'),
                 'bin_files'    : ['_vortexPanel_Cython_cpu.so'],
                 'bin_path'     : os.path.join(kernelsPath,'vortex')},

                {'name'         : 'Source Panel Kernel',
                 'source_type'  : 'Cython',
                 'source_path'  : os.path.join(srcPath,'cython','sourcePanel_Cython'),
                 'bin_files'    : ['_sourcePanel_Cython_cpu.so'],
                 'bin_path'     : os.path.join(kernelsPath,'source')}]   

# panel method functions

print '-----------------------------------------------------------------------\n'
print '-------    P A N E L       M E T H O D      F U N C T I O N S    ------\n'
print '-----------------------------------------------------------------------\n\n'

for kernel in panelKernels:
    # initialize print output
    titleString = kernel['name'] + ' :: ' + kernel['source_type'] + ' code'
    print '== ' + titleString + ' ' + '='*(66-len(titleString))

    # make panelfunctions

    # move to the path of the Makefile
    os.chdir(kernel['source_path'])

    if cleanFlag:
        outputMake = commands.getoutput('make clean')
        print outputMake

    if buildFlag:
        # make file
        outputMake = commands.getoutput('make all')

        # check if binary dir exists
        if not os.path.isdir(kernel['bin_path']):
            os.makedirs(kernel['bin_path'])

        # copy to binary dir
        for binFile in kernel['bin_files']:
            shutil.copy(os.path.join(os.curdir, binFile), kernel['bin_path'])

        print outputMake

    print '======================================================================\n'
