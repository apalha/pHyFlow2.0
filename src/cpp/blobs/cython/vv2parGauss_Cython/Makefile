PYTHON_DIR = /usr/include/python2.7

# define parameters
CC = gcc

all: vv2parGauss_Cython_cpu.so

#--------------------------------------
# 1/r^2 kernel

# ---- Link ---------------------------
vv2parGauss_Cython_cpu.so:  vv2parGauss_Cython_cpu.c
	$(info )
	$(info ************  Compiling  ************)
	$(info )
	$(CC) -fopenmp -shared -pthread -fPIC -fwrapv -O2 -Wall -fno-strict-aliasing -I$(PYTHON_DIR) -o vv2parGauss_Cython_cpu.so vv2parGauss_Cython_cpu.c
 
# ---- gcc compile ------------------
vv2parGauss_Cython_cpu.c:   vv2parGauss_Cython_cpu.pyx
	$(info )
	$(info )
	$(info ************  Generating C file with cython for 1/r^2 kernel ************)
	$(info )
	cython vv2parGauss_Cython_cpu.pyx

#--------------------------------------

clean:
	rm -rf *.o
	rm -rf *.so
	rm -rf vv2parGauss_Cython_cpu.c
