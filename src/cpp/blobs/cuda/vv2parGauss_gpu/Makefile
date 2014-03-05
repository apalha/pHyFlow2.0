PYTHON_DIR = /usr/include/python2.7

# define parameters
CC = gcc
cudaCC = nvcc

all: _vv2parGauss_gpu.so

#-------------------------------------
# Gauss kernel

# ---- Link ---------------------------
_vv2parGauss_gpu.so:  vv2parGauss_gpu.o
	$(info )
	$(info ************  Linking into python module  ************)
	$(info )
	$(CC) -shared -fPIC -O3 -o _vv2parGauss_gpu.so vv2parGauss_gpu.o -L/usr/local/cuda/lib64 -lcudart 
 
# ---- Cuda compile ------------------
vv2parGauss_gpu.o:  vv2parGauss_py_gpu.cu
	$(info )
	$(info )
	$(info ************  COMPILING Gauss kernel ************)
	$(info )
	$(cudaCC) -arch sm_21 -Xcompiler -fPIC -O3 -c vv2parGauss_py_gpu.cu -o vv2parGauss_gpu.o -I$(PYTHON_DIR) -I$(PYTHON_DIR)/numpy

#--------------------------------------

clean:
	rm -rf *.so
	rm -rf *.o
