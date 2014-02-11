PYTHON_DIR = /usr/include/python2.7

all: _interpolation_kernels_1d.so

_interpolation_kernels_1d.so:  _interpolation_kernels_1d.cpp
	$(info )
	$(info ***  Building and linking into python module  ***)
	$(info )
	python setup.py build
	
clean:
	rm -rf *.o
	rm -rf *.so
