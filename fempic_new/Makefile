
# source /ext-home/zl/phd/OP-PIC/scripts/source_oneapi

SRC = .
INC = .
LIB = lib
OBJ = obj
BIN = bin

CPP        = mpicxx
MPICPP	   = mpicxx

CPPFLAGS   = -std=c++11 -Wall -fopenmp -fPIC
NVCCFLAGS  = -m64 -Xptxas -dlcm=ca -Xptxas=-v $(NVCCFLAGS_ADD)

DEBUG ?= 0
ifeq ($(DEBUG), 1)
	CPPFLAGS   += -O0 -g -debug all
	NVCCFLAGS  += -O0 -g -lineinfo
	DEBUG_LOG   = 1
else
	CPPFLAGS   += -O3
	NVCCFLAGS  += -O3
endif

# TODO: Remove; Only for development
CPU_ASSIGN ?= 0
ifeq ($(CPU_ASSIGN), 1)
	NVCCFLAGS  += -DUSE_CPU_ASSIGN
endif

IEEE ?= 0
ifeq ($(IEEE), 1)
	CPPFLAGS   += -fp-model strict -fp-model source -prec-div -prec-sqrt
	NVCCFLAGS  += -fmad=false -ftz=false -prec-div=true -prec-sqrt=true
endif

DEBUG_LOG ?= 0
ifeq ($(DEBUG_LOG), 1)
	CPPFLAGS   += -DDEBUG_LOG
	NVCCFLAGS  += -DDEBUG_LOG
endif

PETSC ?= 1
ifeq ($(PETSC), 1)
	CPPFLAGS   += -DUSE_PETSC
	PETSC_INC  = -I$(PETSC_INSTALL_PATH)/include
	PETSC_LIB  = -L$(PETSC_INSTALL_PATH)/lib -lpetsc
endif

CUDA ?= 1
ifeq ($(CUDA), 1)
	CUDA_INC   = -I$(CUDA_INSTALL_PATH)/include
	CUDA_LIB   = -L$(CUDA_INSTALL_PATH)/lib64 -lcudart -lcuda
endif

ALL_INC    = -I$(OPPIC_PATH)/include $(PETSC_INC) $(CUDA_INC)
ALL_LIBS   = -L$(OPPIC_PATH)/lib $(PETSC_LIB) $(CUDA_LIB)

.PHONY: clean mklib

# ALL_TARGETS = clean fempic_seq fempic_genseq fempic_openmp fempic_genseq_wopet fempic_openmp_wopet fempic_cuda
ALL_TARGETS = clean seq omp cuda

all: $(ALL_TARGETS)

mklib:
	@mkdir -p $(OBJ) $(BIN)

common: mklib 
	     $(CPP) $(CPPFLAGS) -c fempic_op.cpp -o $(OBJ)/fempic_op.o $(ALL_INC)
	     $(CPP) $(CPPFLAGS) -c fempic_ori/meshes.cpp -o $(OBJ)/meshes.o $(ALL_INC)
	     $(CPP) $(CPPFLAGS) -c fempic_ori/particles.cpp -o $(OBJ)/particles.o $(ALL_INC)
		 $(CPP) $(CPPFLAGS) -c fempic_ori/maths.cpp -o $(OBJ)/maths.o $(ALL_INC)
		 $(CPP) $(CPPFLAGS) -c FESolver.cpp -o $(OBJ)/FESolver.o $(ALL_INC)

seq: common
		 $(CPP) $(CPPFLAGS) -c seq/seqkernels.cpp -o $(OBJ)/seqkernels.o $(ALL_INC)
	     $(CPP) $(CPPFLAGS) -o $(BIN)/seq $(OBJ)/fempic_op.o $(OBJ)/meshes.o $(OBJ)/particles.o $(OBJ)/maths.o $(OBJ)/seqkernels.o $(OBJ)/FESolver.o $(ALL_INC) $(ALL_LIBS) -loppic_seq

omp: common
		 $(CPP) $(CPPFLAGS) -c omp/ompkernels.cpp -o $(OBJ)/ompkernels.o $(ALL_INC)
	     $(CPP) $(CPPFLAGS) -o $(BIN)/omp $(OBJ)/fempic_op.o $(OBJ)/meshes.o $(OBJ)/particles.o $(OBJ)/maths.o $(OBJ)/ompkernels.o $(OBJ)/FESolver.o $(ALL_INC) $(ALL_LIBS) -loppic_omp

cuda: common
		 nvcc $(NVCCFLAGS) -c cuda/cudakernels.cu -o $(OBJ)/cudakernels.o $(ALL_INC)
	     $(CPP) $(CPPFLAGS) -o $(BIN)/cuda $(OBJ)/fempic_op.o $(OBJ)/meshes.o $(OBJ)/particles.o $(OBJ)/maths.o $(OBJ)/cudakernels.o $(OBJ)/FESolver.o $(ALL_INC) $(ALL_LIBS) -loppic_cuda


clean: 
	     rm -f *.o *.d *.a
	     rm -f $(OBJ)/*
	     rm -f $(BIN)/*
	     rm -f $(LIB)/*
		 
