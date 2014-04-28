===========
Description
===========

Python Hybrid Flow solver (pHyFlow) is a hybrid flow solver that couples a Navier-Stokes grid solver with a vortex blob solver aimed at solving flows around complex geometries.


============
Installation
============

pHyFlow is based upon several libraries. In order to have a functional version of
pHyFlow follow the following instructions. We assume the user
starts with a fresh installation of ubuntu 12.04LTS. This sequence of instructions
is known to result in a working installation. If you have an already installed ubuntu 12.04LTS
or a more recent one, try following these instructions. If you have any issue please contact us.


1. Install essential packages for compiling programs in ubuntu
	.. code:: python

	sudo apt-get install build-essential

2- Install CUDA
	Go here:
	https://developer.nvidia.com/cuda-downloads?sid=336852

	download Ubuntu 12.04 64bit DEB file

	follow instructions from RPM/DEB installation instructions, basically:

	sudo dpkg -i cuda-repo-<distro>_<version>_<architecture>.deb
	sudo apt-get update
	sudo apt-get install cuda
	
	add the environment variables to cuda in your .bashrc
	
        note that CUDAVERSION is the version of the cuda installed. 5.0, 5.5 or 6.0
        or any other:

	export PATH=/usr/local/cuda-CUDAVERSION/bin:$PATH
	export LD_LIBRARY_PATH=/usr/local/cuda-CUDAVERSION/lib64:$LD_LIBRARY_PATH

	Therefore, for CUDA 6.0:

	export PATH=/usr/local/cuda-6.0/bin:$PATH
	export LD_LIBRARY_PATH=/usr/local/cuda-6.0/lib64:$LD_LIBRARY_PATH

	NOTE: most likely there will be only one CUDA installation in your computer, 
              nevertheless check if in /usr/local there is a folder named cuda with
              the CUDA version you wish to use. If not you will probably have several
              folders, one for each version: cuda-5.0, cuda-5.5, cuda-6.0. In that case,
              you need to make a symbolic link to the version you wish to use. For
              example if you wish to use cuda-6.0 do the following:
	
	      sudo ln -s /usr/local/cuda-6.0 /usr/local/cuda

3- install git, sudo apt-get install git

4- install FEniCS with dorsal
	instructions here: 
		http://fenicsproject.org/download/installation_using_dorsal.html

	or just download:
		https://bitbucket.org/fenics-project/dorsal/get/master.tar.gz

	unpack it, change install directory in dorsal.cfg 
	and use STABLE_BUILD=false and then run:
		./dorsal.sh
		
	follow the instructions, check the software to be installed before
	add 
	
	# FEniCS environment variables
	source /home/bjarnthor/lib/FEniCS/share/fenics/fenics.conf
	
	to .bashrc
5- install boost_python (latest version)
	sudo apt-get install libboost-python-dev

6- install pysparse
	sudo apt-get install python-sparse

7- install setuptools
	sudo apt-get install python-setuptools

8- install pyublas
	Follow the instruction in:
		http://mathema.tician.de/software/pyublas

	or simply download it:
		git clone http://git.tiker.net/trees/pyublas.git

	change directory to the pyublas directory and configure:
		python configure.py

	build pyublas:
		python setup.py build

	install pyublas:
		sudo python setup.py install
	
	copy pyublas include files to /usr/include
		sudo cp -r pyublas/include/pyublas /usr/include

9- install ipython
	sudo apt-get install ipython

10- install scipy
	sudo apt-get install python-scipy

11- install matplotlib
	sudo apt-get install python-matplotlib

12- install spyder
	download the file:
	https://code.google.com/p/spyderlib/downloads/detail?name=spyder-2.2.3.zip

	unzip and run:
		sudo python setup.py install

13- install pHyFlow
	download from bitbucket:
		https://bitbucket.org/gorkiana/phyflow/get/2b2eae42e92a.zip

	extract and change the install.py file to point to the directory desired.

	run python install.py

14- install fenicstools

	download:
		https://github.com/mikaem/fenicstools/archive/master.zip

	extract to a folder named fenictools inside ~/lib and rename it to fenicstools
	
	add folder to PYTHONPATH

15- install pyvtk
	sudo apt-get install python-pyvtk

16- install h5py
	sudo apt-get install python-h5py

17- install mpi4py
	sudo apt-get install python-mpi4py
