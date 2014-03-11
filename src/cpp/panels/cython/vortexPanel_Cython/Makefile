PYTHON_DIR = /usr/include/python2.7

# define parameters
CC = gcc

all: _vortexPanel_Cython_cpu.so

#--------------------------------------
# vortex panels

# ---- Link ---------------------------
_vortexPanel_Cython_cpu.so:  _vortexPanel_Cython_cpu.c
	$(info )
	$(info ************  Compiling  ************)
	$(info )
	$(CC) -fopenmp -shared -pthread -fPIC -fwrapv -O2 -Wall -fno-strict-aliasing -I$(PYTHON_DIR) -o _vortexPanel_Cython_cpu.so _vortexPanel_Cython_cpu.c
 
# ---- gcc compile ------------------
_vortexPanel_Cython_cpu.c:   _vortexPanel_Cython_cpu.pyx
	$(info )
	$(info )
	$(info ************  Generating C file with cython for vortex panels ************)
	$(info )
	cython _vortexPanel_Cython_cpu.pyx

#--------------------------------------

clean:
	rm -rf *.o
	rm -rf *.so
	rm -rf _vortexPanel_Cython_cpu.c
