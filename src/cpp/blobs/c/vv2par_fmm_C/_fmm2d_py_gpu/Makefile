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

CUDA_DEFS = -DCUDASUPPORT

INCLUDE = -I$(CUDA_DIR)/include

LINK = -L$(CUDA_DIR)/lib64


fmm2d_py.o: fmm2d_py.cpp
	$(info )
	$(info = Building :: fmm2d_py.o ====================================)
	$(info )
	$(CC) $(PYTHONFLAGS) $(CXXFLAGS) $(LDFLAGS) -c fmm2d_py.cpp -o fmm2d_py.o $(C_DEFS) $(EXP_DEFS) $(PANEL_DEFS) $(FMM_DEFS) $(CUDA_DEFS) $(INCLUDE)

fmm.o: fmm.cpp
	$(info )
	$(info = Building :: fmm.o =========================================)
	$(info )
	$(CC) $(PYTHONFLAGS) $(CXXFLAGS) $(LDFLAGS) -c fmm.cpp -o fmm.o $(C_DEFS) $(EXP_DEFS) $(PANEL_DEFS) $(FMM_DEFS) $(CUDA_DEFS) $(INCLUDE)

fmmsort.o: fmmsort.cpp
	$(info )
	$(info = Building :: fmmsort.o =====================================)
	$(info )
	$(CC) $(PYTHONFLAGS) $(CXXFLAGS) $(LDFLAGS) -c fmmsort.cpp -o fmmsort.o $(C_DEFS) $(EXP_DEFS) $(PANEL_DEFS) $(FMM_DEFS) $(CUDA_DEFS) $(INCLUDE)

fmmshift.o: fmmshift.cpp
	$(info )
	$(info = Building :: fmmshift.o ====================================)
	$(info )
	$(CC) $(PYTHONFLAGS) $(CXXFLAGS) $(LDFLAGS) -c fmmshift.cpp -o fmmshift.o $(C_DEFS) $(EXP_DEFS) $(PANEL_DEFS) $(FMM_DEFS) $(CUDA_DEFS) $(INCLUDE)

hornershift.o: hornershift.cpp
	$(info )
	$(info = Building :: hornershift.o =================================)
	$(info )
	$(CC) $(PYTHONFLAGS) $(CXXFLAGS) $(LDFLAGS) -c hornershift.cpp -o hornershift.o $(C_DEFS) $(EXP_DEFS) $(PANEL_DEFS) $(FMM_DEFS) $(CUDA_DEFS) $(INCLUDE)

direct.o: direct.cpp
	$(info )
	$(info = Building :: direct.o ======================================)
	$(info )
	$(CC) $(PYTHONFLAGS) $(CXXFLAGS) $(LDFLAGS) -c direct.cpp -o direct.o $(C_DEFS) $(EXP_DEFS) $(PANEL_DEFS) $(FMM_DEFS) $(CUDA_DEFS) $(INCLUDE)

panel.o: panel.cpp
	$(info )
	$(info = Building :: panel.o =======================================)
	$(info )
	$(CC) $(PYTHONFLAGS) $(CXXFLAGS) $(LDFLAGS) -c panel.cpp -o panel.o $(C_DEFS) $(EXP_DEFS) $(PANEL_DEFS) $(FMM_DEFS) $(CUDA_DEFS) $(INCLUDE)

channelpot.o: channelpot.cpp
	$(info )
	$(info = Building :: channelpot.o ==================================)
	$(info )
	$(CC) $(PYTHONFLAGS) $(CXXFLAGS) $(LDFLAGS) -c channelpot.cpp -o channelpot.o $(C_DEFS) $(EXP_DEFS) $(PANEL_DEFS) $(FMM_DEFS) $(CUDA_DEFS) $(INCLUDE)

channelpotpanel.o: channelpotpanel.cpp
	$(info )
	$(info = Building :: channelpotpanel.o =============================)
	$(info )
	$(CC) $(PYTHONFLAGS) $(CXXFLAGS) $(LDFLAGS) -c channelpotpanel.cpp -o channelpotpanel.o $(C_DEFS) $(EXP_DEFS) $(PANEL_DEFS) $(FMM_DEFS) $(CUDA_DEFS) $(INCLUDE)

cudaeval.o: cudaeval.cu
	$(info )
	$(info = Building :: cudaeval.o ====================================)
	$(info )
	$(NVCC) $(LDFLAGSCUDA) $(CXXCUDAFLAGS) -c cudaeval.cu -o cudaeval.o $(C_DEFS) $(EXP_DEFS) $(PANEL_DEFS) $(FMM_DEFS) $(CUDA_DEFS) $(INCLUDE)

cudasort.o: cudasort.cu
	$(info )
	$(info = Building :: cudasort.o ====================================)
	$(info )
	$(NVCC) $(CXXCUDAFLAGS) $(LDFLAGSCUDA) -c cudasort.cu -o cudasort.o $(C_DEFS) $(EXP_DEFS) $(PANEL_DEFS) $(FMM_DEFS) $(CUDA_DEFS) $(INCLUDE)

cudatiming.o: cudatiming.cu
	$(info )
	$(info = Building :: cudatiming.o ==================================)
	$(info )
	$(NVCC) $(CXXCUDAFLAGS) $(LDFLAGSCUDA) -c cudatiming.cu -o cudatiming.o $(C_DEFS) $(EXP_DEFS) $(PANEL_DEFS) $(FMM_DEFS) $(CUDA_DEFS) $(INCLUDE)

_fmm2d_py.so: fmm2d_py.o fmm.o fmmsort.o fmmshift.o hornershift.o direct.o panel.o channelpot.o channelpotpanel.o cudaeval.o cudasort.o cudatiming.o
	$(info )
	$(info = Building and Linking :: _fmm2d_py.so ======================)
	$(info )
	$(CC) -pthread -shared -fPIC -Wl,-O3 -Wl,-Bsymbolic-functions -Wl,-Bsymbolic-functions -Wl,-z,relro fmm2d_py.o fmm.o fmmsort.o fmmshift.o hornershift.o direct.o panel.o channelpot.o channelpotpanel.o -lboost_python cudaeval.o cudasort.o cudatiming.o -o _fmm2d_py.so $(LINK) $(CUDAFLAGS)

all: _fmm2d_py.so

clean:
	rm -rf *.o
	rm -rf *.so

