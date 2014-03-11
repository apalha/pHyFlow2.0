PYTHON_DIR = /usr/include/python2.7

# define parameters
CC = gcc

all: _ww2par_C_cpu.so

#--------------------------------------
# 1/r^2 kernel

# ---- Link ---------------------------
_ww2par_C_cpu.so:  ww2par_C_cpu.o
	$(info )
	$(info ************  Linking into python module  ************)
	$(info )
	$(CC) -fopenmp -shared -fPIC -O2 -o _ww2par_C_cpu.so ww2par_C_cpu.o
 
# ---- gcc compile ------------------
ww2par_C_cpu.o:   ww2par_C_cpu.c
	$(info )
	$(info )
	$(info ************  COMPILING 1/r^2 kernel ************)
	$(info )
	$(CC) -fopenmp -fPIC -O2 -c ww2par_C_cpu.c -o ww2par_C_cpu.o -I$(PYTHON_DIR) -I$(PYTHON_DIR)/numpy

#--------------------------------------

clean:
	rm -rf *.o
	rm -rf *.so
