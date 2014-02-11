PYTHON_DIR = /usr/include/python2.7

# define parameters
CC = gcc
cudaCC = nvcc

all: _vv2par_gpu.so

#-------------------------------------
# 1/r^2 kernel

# ---- Link ---------------------------
_vv2par_gpu.so:  vv2par_gpu.o
	$(info )
	$(info ************  Linking into python module  ************)
	$(info )
	$(CC) -shared -fPIC -O3 -o _vv2par_gpu.so vv2par_gpu.o -L/usr/local/cuda/lib64 -lcudart 
 
# ---- Cuda compile ------------------
vv2par_gpu.o:  vv2par_py_gpu.cu
	$(info )
	$(info )
	$(info ************  COMPILING 1/r^2 kernel ************)
	$(info )
	$(cudaCC) -arch sm_13 -Xcompiler -fPIC -O3 -c vv2par_py_gpu.cu -o vv2par_gpu.o -I$(PYTHON_DIR) -I$(PYTHON_DIR)/numpy

#--------------------------------------

clean:
	rm -rf *.so
	rm -rf *.o
