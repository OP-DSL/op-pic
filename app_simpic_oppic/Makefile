# source /ext-home/zl/phd/OP-PIC/scripts/source_oneapi

SRC = .
INC = .
LIB = lib
OBJ = obj
BIN = bin

CPP        = mpicxx
MPICPP	   = mpicxx

CPPFLAGS   = -std=c++11 -Wall -fopenmp -fPIC
NVCCFLAGS  = -m64 -Xptxas -dlcm=ca -Xptxas=-v

DEBUG ?= 0
ifeq ($(DEBUG), 1)
	CPPFLAGS   += -O0 -g -debug all
	NVCCFLAGS  += -O0 -g -lineinfo
	DEBUG_LOG   = 1
else
	CPPFLAGS   += -O3
	NVCCFLAGS  += -O3
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

TELOS ?= 1
ifeq ($(TELOS), 1)
	NVCCFLAGS  += -gencode arch=compute_70,code=sm_70
	PETSC_PATH = /ext-home/zl/lib_install/petsc-3.18.3-oneapi-2021.3-cuda-rel
else
	NVCCFLAGS  += -gencode arch=compute_60,code=sm_60
	PETSC_PATH = /ext-home/zl/lib_install/petsc-3.18.3-oneapi-2021.3-cuda-kos
endif

PETSC ?= 0
ifeq ($(PETSC), 1)
	CPPFLAGS   += -DUSE_PETSC
	PETSC_INC  = -I$(PETSC_PATH)/include
	PETSC_LIB  = -L$(PETSC_PATH)/lib -lpetsc
else
	PETSC_INC  = 
	PETSC_LIB  = 
endif

CUDA_INSTALL_PATH = /opt/cuda/10.2.89/toolkit
CUDA_INC   = -I$(CUDA_INSTALL_PATH)/include
CUDA_LIB   = -L$(CUDA_INSTALL_PATH)/lib64 -lcudart -lcuda

OPPIC_PATH = /ext-home/zl/phd/OP-PIC/lib_oppic
ALL_INC    = -I$(OPPIC_PATH)/include $(PETSC_INC) $(CUDA_INC)
ALL_LIBS   = -L$(OPPIC_PATH)/lib $(PETSC_LIB) $(CUDA_LIB)

.PHONY: clean mklib

ALL_TARGETS = clean omp genseq seq

all: $(ALL_TARGETS)

mklib:
	@mkdir -p $(OBJ) $(BIN)

common: mklib mklib 
	     $(CPP) $(CPPFLAGS) -c simpic_op.cpp -o $(OBJ)/simpic_op.o $(ALL_INC)
	     $(CPP) $(CPPFLAGS) -c funcs_from_simpic.cpp -o $(OBJ)/funcs_from_simpic.o $(ALL_INC)

genseq: mklib common
		 $(CPP) $(CPPFLAGS) -c seq/seqkernels.cpp -o $(OBJ)/seqkernels.o $(ALL_INC)
	     $(CPP) $(CPPFLAGS) -o $(BIN)/genseq $(OBJ)/simpic_op.o $(OBJ)/funcs_from_simpic.o $(OBJ)/seqkernels.o $(ALL_INC) $(ALL_LIBS) -loppic_seq

omp: mklib common
		 $(CPP) $(CPPFLAGS) -c omp/ompkernels.cpp -o $(OBJ)/ompkernels.o $(ALL_INC)
	     $(CPP) $(CPPFLAGS) -o $(BIN)/omp $(OBJ)/simpic_op.o $(OBJ)/funcs_from_simpic.o $(OBJ)/ompkernels.o $(ALL_INC) $(ALL_LIBS) -loppic_omp

seq: mklib
		$(CPP) $(CPPFLAGS) -c simpic.cpp -o $(OBJ)/simpic.o $(ALL_INC)
		$(CPP) $(CPPFLAGS) -o $(BIN)/seq $(OBJ)/simpic.o $(OBJ)/funcs_from_simpic.o $(ALL_INC) $(ALL_LIBS) -loppic_seq

clean: 
	     rm -f *.o *.d *.a
	     rm -f $(OBJ)/*
	     rm -f $(BIN)/*
	     rm -f $(LIB)/*
