PYTHON_DIR = /usr/include/python2.7
MATLAB_DIR = /usr/local/MATLAB/R2012a
CUDA_DIR = /usr/local/cuda

# define parameters
CC = gcc
cudaCC = nvcc

CXXFLAGS = -O3 -msse3 -fPIC -fno-omit-frame-pointer -pthread

LDFLAGS = -pthread -shared -m64

INCLUDE = -I$(CUDA_DIR)/include

LINK = -L$(CUDA_DIR)/lib64

FMM_DEFS = -DINLINECOMPLEX -DINLINECOMPLEX_DIRECT -DTHETACUTOFFCHECK -DFMM_THETA=0.5 -DFMM_MESHRATIO=2.0 -DCHANNELPOT

PANEL_DEFS = -DPANELSORT -DPANELDUMMYFACTOR=2 -DPANELSHRINKBOX

EXP_DEFS = -DRADIALSHRINK -DRADIALSHRINKLIMIT=512 -DUNROLLHORNER -DUNROLLLEVEL=6 -DSSE2INTRINSIC -DSSE2M2PSSCALING -DSSE2DIRECT -DSSEEVALUNROLL -DSSE2INIT -DFULLCROSSORDER -DOSEENSSEALTERNATIVE

C_DEFS = -DC_CODE

#fmm2dcpu: fmm2d.cpp fmm.cpp fmmsort.cpp fmmshift.cpp hornershift.cpp direct.cpp panel.cpp mexpanel.cpp channelpot.cpp directchannelpot.cpp channelpotpanel.cpp
#	$(CC) $(CXXFLAGS) $(CXXLIBS) $(LDFLAGS) $(INCLUDE) $(LINK) $(FMM_DEFS) $(FMM_DEFS) $(PANEL_DEFS) $(EXP_DEFS) fmm2d.cpp fmm.cpp fmmsort.cpp fmmshift.cpp hornershift.cpp direct.cpp panel.cpp mexpanel.cpp channelpot.cpp directchannelpot.cpp channelpotpanel.cpp

direct.o: direct.cpp
	$(CC) $(C_DEFS) $(CXXFLAGS) $(LDFLAGS) $(INCLUDE) $(LINK) $(FMM_DEFS) $(FMM_DEFS) $(PANEL_DEFS) $(EXP_DEFS) direct.cpp -o direct.o

fmmsort.o: fmmsort.cpp
	$(CC) $(C_DEFS) $(CXXFLAGS) $(LDFLAGS) $(INCLUDE) $(LINK) $(FMM_DEFS) $(FMM_DEFS) $(PANEL_DEFS) $(EXP_DEFS) fmmsort.cpp -o fmmsort.o

fmm.o: fmm.cpp
	$(CC) $(C_DEFS) $(CXXFLAGS) $(LDFLAGS) $(INCLUDE) $(LINK) $(FMM_DEFS) $(FMM_DEFS) $(PANEL_DEFS) $(EXP_DEFS) fmm.cpp -o fmm.o

panel.o: panel.cpp
	$(CC) $(C_DEFS) $(CXXFLAGS) $(LDFLAGS) $(INCLUDE) $(LINK) $(FMM_DEFS) $(FMM_DEFS) $(PANEL_DEFS) $(EXP_DEFS) panel.cpp -o panel.o

fmmshift.o: fmmshift.cpp
	$(CC) $(C_DEFS) $(CXXFLAGS) $(LDFLAGS) $(INCLUDE) $(LINK) $(FMM_DEFS) $(FMM_DEFS) $(PANEL_DEFS) $(EXP_DEFS) fmmshift.cpp -o fmmshift.o

hornershift.o: hornershift.cpp
	$(CC) $(C_DEFS) $(CXXFLAGS) $(LDFLAGS) $(INCLUDE) $(LINK) $(FMM_DEFS) $(FMM_DEFS) $(PANEL_DEFS) $(EXP_DEFS) hornershift.cpp -o hornershift.o

channelpot.o: channelpot.cpp
	$(CC) $(C_DEFS) $(CXXFLAGS) $(LDFLAGS) $(INCLUDE) $(LINK) $(FMM_DEFS) $(FMM_DEFS) $(PANEL_DEFS) $(EXP_DEFS) channelpot.cpp -o channelpot.o

channelpotpanel.o: channelpotpanel.cpp
	$(CC) $(C_DEFS) $(CXXFLAGS) $(LDFLAGS) $(INCLUDE) $(LINK) $(FMM_DEFS) $(FMM_DEFS) $(PANEL_DEFS) $(EXP_DEFS) channelpotpanel.cpp -o channelpotpanel.o

directchannelpot.o: directchannelpot.cpp
	$(CC) $(C_DEFS) $(CXXFLAGS) $(LDFLAGS) $(INCLUDE) $(LINK) $(FMM_DEFS) $(FMM_DEFS) $(PANEL_DEFS) $(EXP_DEFS) directchannelpot.cpp -o directchannelpot.o

all: direct.o fmmsort.o fmm.o panel.o fmmshift.o hornershift.o channelpot.o channelpotpanel.o directchannelpot.o

clean:
	rm -rf *.o
