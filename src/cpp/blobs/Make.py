"""Python Makefile to automate the build and installation of vortex C Python modules
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

import commands
import os
import sys
import shutil


kernelsPath = os.path.join(externCodeBlobsPath,'kernels')
remeshPath = os.path.join(externCodeBlobsPath,'remesh')
srcPath = os.path.abspath(os.path.curdir)

# check if libPath already exists if not exit
if not os.path.isdir(libRoot):
    raise BaseException('library path does not exist')

# induced velocity kernels
vkernels = [{'name':'Cutoff Kernel', 'source_type':'C', 'source_path':os.path.join(srcPath,'c','vv2par_C'),'bin_files':['_vv2par_C_cpu.so'],'bin_path':os.path.join(kernelsPath,'cutoff'),'make_options':''},\
            {'name':'Cutoff Kernel', 'source_type':'Cython', 'source_path':os.path.join(srcPath,'cython','vv2par_Cython'),'bin_files':['vv2par_Cython_cpu.so'],'bin_path':os.path.join(kernelsPath,'cutoff'),'make_options':''},\
            {'name':'Cutoff Kernel', 'source_type':'GPU', 'source_path':os.path.join(srcPath,'cuda','vv2par_gpu'),'bin_files':['_vv2par_gpu.so'],'bin_path':os.path.join(kernelsPath,'cutoff'),'make_options':''},\
            {'name':'Gauss Kernel', 'source_type':'C', 'source_path':os.path.join(srcPath,'c','vv2parGauss_C'),'bin_files':['_vv2parGauss_C_cpu.so'],'bin_path':os.path.join(kernelsPath,'gauss'),'make_options':''},\
            {'name':'Gauss Kernel', 'source_type':'Cython', 'source_path':os.path.join(srcPath,'cython','vv2parGauss_Cython'),'bin_files':['vv2parGauss_Cython_cpu.so'],'bin_path':os.path.join(kernelsPath,'gauss'),'make_options':''},\
            {'name':'Gauss Kernel', 'source_type':'GPU', 'source_path':os.path.join(srcPath,'cuda','vv2parGauss_gpu'),'bin_files':['_vv2parGauss_gpu.so'],'bin_path':os.path.join(kernelsPath,'gauss'),'make_options':''},\
            {'name':'Gauss Kernel (FMM/CPU)', 'source_type':'C/C++', 'source_path':os.path.join(srcPath,'c++','fmm2d_py_cpu_gpu'),'bin_files':['_fmm2d_py_cpu.so'],'bin_path':os.path.join(kernelsPath,'gauss'),'make_options':''},\
            {'name':'Gauss Kernel (FMM/GPU)', 'source_type':'C/C++/CUDA', 'source_path':os.path.join(srcPath,'c++','fmm2d_py_cpu_gpu'),'bin_files':['_fmm2d_py_gpu.so'],'bin_path':os.path.join(kernelsPath,'gauss'),'make_options':'CUDA_DEFS=-DCUDASUPPORT'}

]

# induced vorticity kernels
wkernels = [{'name':'Cutoff Kernel', 'source_type':'C', 'source_path':os.path.join(srcPath,'c','ww2par_C'),'bin_files':['_ww2par_C_cpu.so'],'bin_path':os.path.join(kernelsPath,'cutoff')},\
            {'name':'Cutoff Kernel', 'source_type':'GPU', 'source_path':os.path.join(srcPath,'cuda','ww2par_gpu'),'bin_files':['_ww2par_gpu.so'],'bin_path':os.path.join(kernelsPath,'cutoff')},\
            {'name':'Gauss Kernel', 'source_type':'C', 'source_path':os.path.join(srcPath,'c','ww2parGauss_C'),'bin_files':['_ww2parGauss_C_cpu.so'],'bin_path':os.path.join(kernelsPath,'gauss')},\
            {'name':'Gauss Kernel', 'source_type':'GPU', 'source_path':os.path.join(srcPath,'cuda','ww2parGauss_gpu'),'bin_files':['_ww2parGauss_gpu.so'],'bin_path':os.path.join(kernelsPath,'gauss')}]

# interpolation kernels  (used for remeshing vortex blobs)
interpolationkernels = [{'name':'Interpolation kernels', 'source_type':'C++', 'source_path':os.path.join(srcPath,'c++','interpolation_kernels_1d'),'bin_files':['_interpolation_kernels_1d.so'],'bin_path':remeshPath}]


# update the velocity vortex functions

print '-----------------------------------------------------------------------\n'
print '-------       V E L O C I T Y               K E R N E L S       -------\n'
print '-----------------------------------------------------------------------\n\n'

for vortexKernel in vkernels:
    # initialize print output
    titleString = vortexKernel['name'] + ' :: ' + vortexKernel['source_type'] + ' code'
    print '== ' + titleString + ' ' + '='*(66-len(titleString))
    
    # move to the path of the Makefile
    os.chdir(vortexKernel['source_path'])
        
    if cleanFlag:
        outputMake = commands.getoutput('make clean')
        print outputMake
    
    # make kernels
    
    if buildFlag:    
        # make files
        outputMake = commands.getoutput('make all ' + vortexKernel['make_options'])
        
        # check if binary dir exists        
        if not os.path.isdir(vortexKernel['bin_path']):
            os.makedirs(vortexKernel['bin_path'])
        
        # copy to binary dir
        for binFile in vortexKernel['bin_files']:
            shutil.copy(os.path.join(os.curdir,binFile),vortexKernel['bin_path'])
        
        print outputMake
        
    print '======================================================================\n'



# update the vorticity vortex functions

print '-----------------------------------------------------------------------\n'
print '-------       V O R T I C I T Y             K E R N E L S       -------\n'
print '-----------------------------------------------------------------------\n\n'

for vortexKernel in wkernels:
    # initialize print output
    titleString = vortexKernel['name'] + ' :: ' + vortexKernel['source_type'] + ' code'
    print '== ' + titleString + ' ' + '='*(66-len(titleString))
    
    # move to the path of the Makefile
    os.chdir(vortexKernel['source_path'])
        
    if cleanFlag:
        outputMake = commands.getoutput('make clean')
        print outputMake
    
    # make kernels
    
    if buildFlag:    
        # make files
        outputMake = commands.getoutput('make all')
        
        # check if binary dir exists        
        if not os.path.isdir(vortexKernel['bin_path']):
            os.makedirs(vortexKernel['bin_path'])
        
        # copy to binary dir
        for binFile in vortexKernel['bin_files']:
            shutil.copy(os.path.join(os.curdir,binFile),vortexKernel['bin_path'])
        
        print outputMake
        
    print '======================================================================\n'



# update the resmeshing interpolation kernels

print '-----------------------------------------------------------------------\n'
print '-------       R E M E S H I N G             K E R N E L S       -------\n'
print '-----------------------------------------------------------------------\n\n'

for interpKernel in interpolationkernels:
    # initialize print output
    titleString = interpKernel['name'] + ' :: ' + interpKernel['source_type'] + ' code'
    print '== ' + titleString + ' ' + '='*(66-len(titleString))
    
    # move to the path of the Makefile
    os.chdir(interpKernel['source_path'])
        
    if cleanFlag:
        outputMake = commands.getoutput('make clean')
        print outputMake
    
    # make kernels
    
    if buildFlag:    
        # make files
        outputMake = commands.getoutput('make all')
        
        # check if binary dir exists        
        if not os.path.isdir(interpKernel['bin_path']):
            os.makedirs(interpKernel['bin_path'])
        
        # copy to binary dir
        for binFile in interpKernel['bin_files']:
            shutil.copy(os.path.join(os.curdir,binFile),interpKernel['bin_path'])
        
        print outputMake
        
    print '======================================================================\n'
