PYTHON_DIR = /usr/include/python2.7

# define parameters
CC = gcc

all: _ww2parGauss_C_cpu.so

#--------------------------------------
# Gaussian kernel

# ---- Link ---------------------------
_ww2parGauss_C_cpu.so:  ww2parGauss_C_cpu.o
	$(info )
	$(info ************  Linking into python module  ************)
	$(info )
	$(CC) -fopenmp -shared -fPIC -O3 -o _ww2parGauss_C_cpu.so ww2parGauss_C_cpu.o
 
# ---- gcc compile ------------------
ww2parGauss_C_cpu.o:   ww2parGauss_C_cpu.c
	$(info )
	$(info )
	$(info ************  COMPILING Gaussian kernel  ************)
	$(info )
	$(CC) -fopenmp -fPIC -O3 -c ww2parGauss_C_cpu.c -o ww2parGauss_C_cpu.o -I$(PYTHON_DIR) -I$(PYTHON_DIR)/numpy

#------------------------------------

clean:
	rm -rf *.o
	rm -rf *.so
