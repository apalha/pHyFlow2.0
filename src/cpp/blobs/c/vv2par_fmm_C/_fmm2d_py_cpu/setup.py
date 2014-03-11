#!/usr/bin/env python

from distutils.core import setup
from distutils.extension import Extension
import os
import shutil

os.environ["CC"] = "g++"


INCLUDE = "-I/usr/local/cuda/include"


setup(name="_fmm2d_py",
    ext_modules=[
        Extension("_fmm2d_py", ["fmm2d_py.cpp","fmm.cpp", "fmmsort.cpp","fmmshift.cpp","hornershift.cpp","direct.cpp","panel.cpp", "channelpot.cpp","channelpotpanel.cpp"],\
        libraries = ["boost_python"],\
        extra_compile_args = ["-O3", "-msse3", "-fPIC", "-fno-omit-frame-pointer", "-pthread", "-shared", "-m64", "-DC_CODE", "-DRADIALSHRINK", "-DRADIALSHRINKLIMIT=512", "-DUNROLLHORNER", "-DUNROLLLEVEL=6", "-DSSE2INTRINSIC",\
                              "-DSSE2M2PSSCALING", "-DSSE2DIRECT", "-DSSEEVALUNROLL", "-DSSE2INIT", "-DFULLCROSSORDER", "-DOSEENSSEALTERNATIVE", "-DPANELSORT",\
                              "-DPANELDUMMYFACTOR=2", "-DPANELSHRINKBOX", "-DINLINECOMPLEX", "-DINLINECOMPLEX_DIRECT", "-DTHETACUTOFFCHECK", "-DFMM_THETA=0.5", "-DFMM_MESHRATIO=2.0", "-DCHANNELPOT", "-I/usr/local/cuda/include"],\
        extra_link_args=["-L/usr/local/cuda/lib64"])])

print '************  Copying compiled file  ************'

# copy compiled extension to root directory
list_dir = os.listdir(os.path.join(os.curdir,'build')) # get the list of build directories
build_dir = [myDir for myDir in list_dir if "lib." in myDir] # find the ones that start with lib.

if len(build_dir)>1:
    raise NameError('Too many build directories! Clean build directory first!')

# copy the file
shutil.copy(os.path.join(os.curdir,'build',build_dir[0],'_fmm2d_py.so'),os.curdir)


