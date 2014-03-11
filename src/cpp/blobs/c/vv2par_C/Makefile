PYTHON_DIR = /usr/include/python2.7

# define parameters
CC = gcc

all: _vv2par_C_cpu.so

#--------------------------------------
# 1/r^2 kernel

# ---- Link ---------------------------
_vv2par_C_cpu.so:  vv2par_C_cpu.o
	$(info )
	$(info ************  Linking into python module  ************)
	$(info )
	$(CC) -fopenmp -shared -fPIC -O2 -o _vv2par_C_cpu.so vv2par_C_cpu.o
 
# ---- gcc compile ------------------
vv2par_C_cpu.o:   vv2par_C_cpu.c
	$(info )
	$(info )
	$(info ************  COMPILING 1/r^2 kernel ************)
	$(info )
	$(CC) -fopenmp -fPIC -O2 -c vv2par_C_cpu.c -o vv2par_C_cpu.o -I$(PYTHON_DIR) -I$(PYTHON_DIR)/numpy

#--------------------------------------

clean:
	rm -rf *.o
	rm -rf *.so
