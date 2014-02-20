PYTHON_DIR = /usr/include/python2.7

# define parameters
CC = gcc

all: _sourcePanel_Cython_cpu.so

#--------------------------------------
# vortex panels

# ---- Link ---------------------------
_sourcePanel_Cython_cpu.so:  _sourcePanel_Cython_cpu.c
	$(info )
	$(info ************  Compiling  ************)
	$(info )
	$(CC) -fopenmp -shared -pthread -fPIC -fwrapv -O2 -Wall -fno-strict-aliasing -I$(PYTHON_DIR) -o _sourcePanel_Cython_cpu.so _sourcePanel_Cython_cpu.c
 
# ---- gcc compile ------------------
_sourcePanel_Cython_cpu.c:   _sourcePanel_Cython_cpu.pyx
	$(info )
	$(info )
	$(info ************  Generating C file with cython for source panels ************)
	$(info )
	cython _sourcePanel_Cython_cpu.pyx

#--------------------------------------

clean:
	rm -rf *.o
	rm -rf *.so
	rm -rf _sourcePanel_Cython_cpu.c
