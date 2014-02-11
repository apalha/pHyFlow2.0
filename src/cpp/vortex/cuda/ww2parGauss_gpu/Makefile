PYTHON_DIR = /usr/include/python2.7

# define parameters
CC = gcc
cudaCC = nvcc

all: _ww2parGauss_gpu.so

#-------------------------------------
# Gauss kernel

# ---- Link ---------------------------
_ww2parGauss_gpu.so:  ww2parGauss_gpu.o
	$(info )
	$(info ************  Linking into python module  ************)
	$(info )
	$(CC) -shared -fPIC -O3 -o _ww2parGauss_gpu.so ww2parGauss_gpu.o -L/usr/local/cuda/lib64 -lcudart 
 
# ---- Cuda compile ------------------
ww2parGauss_gpu.o:  ww2parGauss_py_gpu.cu
	$(info )
	$(info )
	$(info ************  COMPILING Gauss kernel ************)
	$(info )
	$(cudaCC) -arch sm_21 -Xcompiler -fPIC -O3 -c ww2parGauss_py_gpu.cu -o ww2parGauss_gpu.o -I$(PYTHON_DIR) -I$(PYTHON_DIR)/numpy

#--------------------------------------

clean:
	rm -rf *.so
	rm -rf *.o
