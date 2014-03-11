PYTHON_DIR = /usr/include/python2.7

# define parameters
CC = gcc
cudaCC = nvcc

all: _ww2par_gpu.so

#-------------------------------------
# 1/r^2 kernel

# ---- Link ---------------------------
_ww2par_gpu.so:  ww2par_gpu.o
	$(info )
	$(info ************  Linking into python module  ************)
	$(info )
	$(CC) -shared -fPIC -O3 -o _ww2par_gpu.so ww2par_gpu.o -L/usr/local/cuda/lib64 -lcudart 
 
# ---- Cuda compile ------------------
ww2par_gpu.o:  ww2par_py_gpu.cu
	$(info )
	$(info )
	$(info ************  COMPILING 1/r^2 kernel ************)
	$(info )
	$(cudaCC) -arch sm_13 -Xcompiler -fPIC -O3 -c ww2par_py_gpu.cu -o ww2par_gpu.o -I$(PYTHON_DIR) -I$(PYTHON_DIR)/numpy

#--------------------------------------

clean:
	rm -rf *.so
	rm -rf *.o
