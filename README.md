# pHyFlow

## Description

Python Hybrid Flow solver (pHyFlow) is a hybrid flow solver that couples a Navier-Stokes grid solver with a vortex blob solver aimed at solving flows around complex geometries.


## Installation

pHyFlow is based upon several libraries. In order to have a functional version of
pHyFlow follow the following instructions. We assume the user
starts with a fresh installation of ubuntu 12.04LTS. This sequence of instructions
is known to result in a working installation. If you have an already installed ubuntu 12.04LTS
or a more recent one, try following these instructions. If you have any issue please contact us.


### 1. Install essential packages for compiling programs in ubuntu
In the terminal line type:

	sudo apt-get install build-essential

### 2. Install CUDA
Go to this link: [https://developer.nvidia.com/cuda-downloads?sid=336852](https://developer.nvidia.com/cuda-downloads?sid=336852) and
download the Ubuntu 12.04 64bit DEB file or the one corresponding to your Ubuntu version.
Follow the DEB installation instructions, basically:

	sudo dpkg -i cuda-repo-<distro>_<version>_<architecture>.deb
	sudo apt-get update
	sudo apt-get install cuda
	
Add the environment variables to cuda in your `.bashrc`:
	
	export PATH=/usr/local/cuda-CUDAVERSION/bin:$PATH
	export LD_LIBRARY_PATH=/usr/local/cuda-CUDAVERSION/lib64:$LD_LIBRARY_PATH

Note that CUDAVERSION is the version of the cuda installed. 5.0, 5.5 or 6.0 or any other.
Therefore, for CUDA 6.0:

	export PATH=/usr/local/cuda-6.0/bin:$PATH
	export LD_LIBRARY_PATH=/usr/local/cuda-6.0/lib64:$LD_LIBRARY_PATH

NOTE: most likely there will be only one CUDA installation in your computer, 
nevertheless check if in /usr/local there is a folder named cuda with
the CUDA version you wish to use. If not you will probably have several
folders, one for each version: `cuda-5.0`, `cuda-5.5`, `cuda-6.0`. In that case,
you need to make a symbolic link to the version you wish to use. For
example if you wish to use `cuda-6.0` do the following:
	
      sudo ln -s /usr/local/cuda-6.0 /usr/local/cuda

### 3. Install git
In the terminal line type:
	
	sudo apt-get install git

### 4. Install FEniCS using dorsal
You have two options: you can follow dorsal's instructions in [http://fenicsproject.org/download/installation_using_dorsal.html](http://fenicsproject.org/download/installation_using_dorsal.html)
or just download [https://bitbucket.org/fenics-project/dorsal/get/master.tar.gz](https://bitbucket.org/fenics-project/dorsal/get/master.tar.gz)
unpack it, change the install directory in `dorsal.cfg` to the one you wish, for example `install_dir_FEniCS`,
and use the option

	STABLE_BUILD=false

then run:

	./dorsal.sh
		
*IMPORTANT*: follow the instructions, check the software to be installed before continuing with dorsal's installation.
	
Add the environment variables to FEniCS in your `.bashrc`:
	
	# FEniCS environment variables
	source install_dir_FEniCS/FEniCS/share/fenics/fenics.conf

### 5. Install boost python (latest version)
In the terminal line type:

	sudo apt-get install libboost-python-dev

### 6. Install pysparse
In the terminal line type:

	sudo apt-get install python-sparse

### 7. Install setuptools
In the terminal line type:	
	
	sudo apt-get install python-setuptools

### 8. Install pyublas
You can follow the instructions in [http://mathema.tician.de/software/pyublas](http://mathema.tician.de/software/pyublas)
or simply download it:

	git clone http://git.tiker.net/trees/pyublas.git

then change directory to the pyublas directory and configure:

	python configure.py

build pyublas:

	python setup.py build

install pyublas:

	sudo python setup.py install
	
copy pyublas include files to `/usr/include`

	sudo cp -r pyublas/include/pyublas /usr/include

### 9. Install ipython
In the terminal line type:
	
	sudo apt-get install ipython

### 10. Install scipy
In the terminal line type:
	
	sudo apt-get install python-scipy

### 11. Install matplotlib
In the terminal line type:

	sudo apt-get install python-matplotlib

### 12. Install spyder (optional)
Download the latest version from [https://bitbucket.org/spyder-ide/spyderlib/downloads](https://bitbucket.org/spyder-ide/spyderlib/downloads) unzip it and run:

	sudo python setup.py install

As an alternative you can follow the instructions here: [https://pythonhosted.org/spyder/installation.html](https://pythonhosted.org/spyder/installation.html).

### 13. Install pHyFlow
Create a directory for your pHyFlow download folder. Move into that folder and download the latest version of pHyFlow using git:

	git clone https://gorkiana@bitbucket.org/gorkiana/phyflow2.0.git

Change in the `install.py` file the line

	libParent = '/media/DATAPART1/programs/lib/python' # root install directory

into whatever directory you wish to use to install pHyFlow, for example `/home/your_username/intall_pHyFlow_parent_dir/`. Then type on the terminal line:

	run python install.py

Add to your `.bashrc` the pHyFlow install directory in order to be able to import pHyFlow into python:
	
	# add pHyFlow to the python path                                                                                      
	export PYTHONPATH=$PYTHONPATH:/home/your_username/intall_pHyFlow_parent_dir
	

### 14. Install fenicstools
Download the lastest release of fenicstools from [https://github.com/mikaem/fenicstools/releases](https://github.com/mikaem/fenicstools/releases)
Extract it into `/home/your_username/intall_pHyFlow_parent_dir/`, as used in (13) and rename the extract folder to fenicstools. Now you do not need
to add this folder to `PYTHONPATH` again since it was already done in (13).

### 15. Install pyvtk
In the terminal line type:

	sudo apt-get install python-pyvtk

### 16. Install h5py
In the terminal line type:
	
	sudo apt-get install python-h5py

### 17. Install mpi4py
In the terminal line type:

	sudo apt-get install python-mpi4py

