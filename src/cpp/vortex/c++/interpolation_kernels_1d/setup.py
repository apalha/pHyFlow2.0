#!/usr/bin/env python
 
from distutils.core import setup
from distutils.extension import Extension
import os
import shutil

setup(name="_interpolation_kernels_1d",
    ext_modules=[
        Extension("_interpolation_kernels_1d", ["_interpolation_kernels_1d.cpp"],
        libraries = ["boost_python"],
        extra_compile_args = ["-fopenmp"],
        extra_link_args=['-lgomp'])
    ])

print '************  Copying compiled file  ************'

# copy compiled extension to root directory
list_dir = os.listdir(os.path.join(os.curdir,'build')) # get the list of build directories
build_dir = [myDir for myDir in list_dir if "lib." in myDir] # find the ones that start with lib.

if len(build_dir)>1:
    raise NameError('Too many build directories! Clean build directory first!')

# copy the file
shutil.copy(os.path.join(os.curdir,'build',build_dir[0],'_interpolation_kernels_1d.so'),os.curdir)


