
ifeq ($(strip $(CUDA_DEFS)),-DCUDASUPPORT)
# GPU definitions

# define the directories where libraries lie
PYTHON_DIR = /usr/include/python2.7
CUDA_DIR = /usr/local/cuda

# define compile parameters

# compilers
CC = g++
NVCC = nvcc

# python flags required to generate the python module
PYTHONFLAGS = -fno-strict-aliasing -DNDEBUG -g -fwrapv -I/usr/include/python2.7

# compiler optimization flags
CXXFLAGS = -w -O3 -msse3 -fno-omit-frame-pointer
LDFLAGS = -pthread -fPIC -shared -m64
LDFLAGSCUDA = -arch sm_21 -Xcompiler -fPIC
CXXCUDAFLAGS = -w -O3
CUDAFLAGS = -lcudart
FMM_DEFS = -DINLINECOMPLEX -DINLINECOMPLEX_DIRECT -DTHETACUTOFFCHECK -DFMM_THETA=0.5 -DFMM_MESHRATIO=2.0 -DCHANNELPOT
PANEL_DEFS = -DPANELSORT -DPANELDUMMYFACTOR=2 -DPANELSHRINKBOX
EXP_DEFS = -DRADIALSHRINK -DRADIALSHRINKLIMIT=512 -DUNROLLHORNER -DUNROLLLEVEL=6 -DSSE2INTRINSIC -DSSE2M2PSSCALING -DSSE2DIRECT -DSSEEVALUNROLL -DSSE2INIT -DFULLCROSSORDER -DOSEENSSEALTERNATIVE
C_DEFS = -DC_CODE
CUDA_DEFS = -DCUDASUPPORT -DMULTIPOLEEVALCUDA -DMULTIPOLEINITCUDA -DMULTIPOLESHIFTCUDA -DCUDASORT
INCLUDE = -I$(CUDA_DIR)/include
LINK = -L$(CUDA_DIR)/lib64

# compilation of the python module (interface to C/C++/CUDA)

# CPP files compilation
%_gpu.o: %.cpp
	$(info )
	$(info = Building :: $*_gpu.o ====================================)
	$(info )
	$(CC) $(PYTHONFLAGS) $(CXXFLAGS) $(LDFLAGS) -c $*.cpp -o $*_gpu.o $(C_DEFS) $(EXP_DEFS) $(PANEL_DEFS) $(FMM_DEFS) $(CUDA_DEFS) $(INCLUDE)

# CUDA files compilation
%_gpu.o: %.cu
	$(info )
	$(info = Building :: $*.o ====================================)
	$(info )
	$(NVCC) $(LDFLAGSCUDA) $(CXXCUDAFLAGS) -c $*.cu -o $*_gpu.o $(C_DEFS) $(EXP_DEFS) $(PANEL_DEFS) $(FMM_DEFS) $(CUDA_DEFS) $(INCLUDE)

# linking of the python module
_fmm2d_py_gpu.so: fmm2d_py_gpu.o fmm_gpu.o fmmsort_gpu.o fmmshift_gpu.o hornershift_gpu.o direct_gpu.o panel_gpu.o channelpot_gpu.o channelpotpanel_gpu.o cudaeval_gpu.o cudasort_gpu.o cudatiming_gpu.o
	$(info )
	$(info = Building and Linking :: _fmm2d_py_gpu.so ======================)
	$(info )
	$(CC) -pthread -shared -fPIC -Wl,-O3 -Wl,-Bsymbolic-functions -Wl,-Bsymbolic-functions -Wl,-z,relro fmm2d_py_gpu.o fmm_gpu.o fmmsort_gpu.o fmmshift_gpu.o hornershift_gpu.o direct_gpu.o panel_gpu.o channelpot_gpu.o channelpotpanel_gpu.o -lboost_python cudaeval_gpu.o cudasort_gpu.o cudatiming_gpu.o -o _fmm2d_py_gpu.so $(LINK) $(CUDAFLAGS)

all: _fmm2d_py_gpu.so

else
# CPU definitions

PYTHON_DIR = /usr/include/python2.7

# define compile parameters
CC = g++

PYTHONFLAGS = -fno-strict-aliasing -DNDEBUG -g -fwrapv -I/usr/include/python2.7

CXXFLAGS = -w -O3 -msse3 -fPIC -fno-omit-frame-pointer

LDFLAGS = -pthread -shared -m64

FMM_DEFS = -DINLINECOMPLEX -DINLINECOMPLEX_DIRECT -DTHETACUTOFFCHECK -DFMM_THETA=0.5 -DFMM_MESHRATIO=2.0 -DCHANNELPOT

PANEL_DEFS = -DPANELSORT -DPANELDUMMYFACTOR=2 -DPANELSHRINKBOX

EXP_DEFS = -DRADIALSHRINK -DRADIALSHRINKLIMIT=512 -DUNROLLHORNER -DUNROLLLEVEL=6 -DSSE2INTRINSIC -DSSE2M2PSSCALING -DSSE2DIRECT -DSSEEVALUNROLL -DSSE2INIT -DFULLCROSSORDER -DOSEENSSEALTERNATIVE

C_DEFS = -DC_CODE

# CPP files compilation
%_cpu.o: %.cpp
	$(info )
	$(info = Building :: $*_cpu.o ====================================)
	$(info )
	$(CC) $(PYTHONFLAGS) $(CXXFLAGS) $(LDFLAGS) -c $*.cpp -o $*_cpu.o $(C_DEFS) $(EXP_DEFS) $(PANEL_DEFS) $(FMM_DEFS)

_fmm2d_py_cpu.so: fmm2d_py_cpu.o fmm_cpu.o fmmsort_cpu.o fmmshift_cpu.o hornershift_cpu.o direct_cpu.o panel_cpu.o channelpot_cpu.o channelpotpanel_cpu.o
	$(info )
	$(info = Building and Linking :: _fmm2d_py_cpu.so ======================)
	$(info )
	$(CC) -pthread -shared -Wl,-O3 -Wl,-Bsymbolic-functions -Wl,-Bsymbolic-functions -Wl,-z,relro fmm2d_py_cpu.o fmm_cpu.o fmmsort_cpu.o fmmshift_cpu.o hornershift_cpu.o direct_cpu.o panel_cpu.o channelpot_cpu.o channelpotpanel_cpu.o -lboost_python -o _fmm2d_py_cpu.so

all: _fmm2d_py_cpu.so

endif


clean:
	rm -rf *.o
	rm -rf *.so

clean_cpu:
	rm -rf *_cpu.o
	rm -rf *_cpu.so

clean_gpu:
	rm -rf *_gpu.o
	rm -rf *_gpu.so

