PYTHON_DIR = /usr/include/python2.7

# define parameters
CC = gcc

all: _vortexPanel_Cython_cpu_normal.so

#--------------------------------------
# vortex panels

# ---- Link ---------------------------
_vortexPanel_Cython_cpu_normal.so:  _vortexPanel_Cython_cpu_normal.c
	$(info )
	$(info ************  Compiling  ************)
	$(info )
	$(CC) -fopenmp -shared -pthread -fPIC -fwrapv -O2 -Wall -fno-strict-aliasing -I$(PYTHON_DIR) -o _vortexPanel_Cython_cpu_normal.so _vortexPanel_Cython_cpu_normal.c
 
# ---- gcc compile ------------------
_vortexPanel_Cython_cpu_normal.c:   _vortexPanel_Cython_cpu_normal.pyx
	$(info )
	$(info )
	$(info ************  Generating C file with cython for vortex panels ************)
	$(info )
	cython _vortexPanel_Cython_cpu_normal.pyx

#--------------------------------------

clean:
	rm -rf *.o
	rm -rf *.so
	rm -rf _vortexPanel_Cython_cpu_normal.c
