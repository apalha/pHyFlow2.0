PYTHON_DIR = /usr/include/python2.7

# define parameters
CC = gcc

all: _vv2parGauss_C_cpu.so

#--------------------------------------
# Gaussian kernel

# ---- Link ---------------------------
_vv2parGauss_C_cpu.so:  vv2parGauss_C_cpu.o
	$(info )
	$(info ************  Linking into python module  ************)
	$(info )
	$(CC) -fopenmp -shared -fPIC -O3 -o _vv2parGauss_C_cpu.so vv2parGauss_C_cpu.o
 
# ---- gcc compile ------------------
vv2parGauss_C_cpu.o:   vv2parGauss_C_cpu.c
	$(info )
	$(info )
	$(info ************  COMPILING Gaussian kernel  ************)
	$(info )
	$(CC) -fopenmp -fPIC -O3 -c vv2parGauss_C_cpu.c -o vv2parGauss_C_cpu.o -I$(PYTHON_DIR) -I$(PYTHON_DIR)/numpy

#------------------------------------

clean:
	rm -rf *.o
	rm -rf *.so
