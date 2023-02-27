
# module load gnu-6.3.0/gnu-6.3.0
# module load gnu-6.3.0/ompi-1.10.6
# module load cuda/toolkit-10.2.89

# module load intel-oneapi-2021.3/intel-classic-2021.3 intel-oneapi-2021.3/impi-classic-2021.3
# module load cuda/toolkit-10.2.89

# export LD_LIBRARY_PATH=/ext-home/zl/phd/OP-PIC/app_fempic_oppic/lib:/ext-home/zl/lib_install/petsc-3.18.2-oneapi-2021.3/lib:$LD_LIBRARY_PATH
# export LD_LIBRARY_PATH=/ext-home/zl/phd/OP-PIC/app_fempic_oppic/lib:/ext-home/zl/lib_install/petsc-3.18.2-oneapi-2021.3-release/lib:$LD_LIBRARY_PATH

SRC = src
INC = include
LIB = lib
OBJ = obj

CPP        = mpicxx
MPICPP	   = mpicxx

CPPFLAGS   = -std=c++11 -Wall -fopenmp -fPIC
NVCCFLAGS  = -m64 -Xptxas -dlcm=ca -Xptxas=-v

DEBUG ?= 0
ifeq ($(DEBUG), 1)
	CPPFLAGS   += -O0 -g
	NVCCFLAGS  += -O0 -g
else
	CPPFLAGS   += -O3
	NVCCFLAGS  += -O3
endif

TELOS ?= 1
ifeq ($(TELOS), 1)
	NVCCFLAGS  += -gencode arch=compute_70,code=sm_70
	PETSC_PATH = /ext-home/zl/lib_install/petsc-3.18.3-oneapi-2021.3-cuda-rel
else
	NVCCFLAGS  += -gencode arch=compute_60,code=sm_60
	PETSC_PATH = /ext-home/zl/lib_install/petsc-3.18.3-oneapi-2021.3-cuda-kos
endif

PETSC ?= 1
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

ALL_TARGETS = clean oppic_seq oppic_omp oppic_cuda

all: $(ALL_TARGETS)

mklib:
	@mkdir -p $(LIB) $(OBJ)


oppic_core: mklib
	     $(CPP) $(CPPFLAGS) -c $(OPPIC_PATH)/$(SRC)/core/oppic_lib_core.cpp -o $(OBJ)/oppic_lib_core.o $(ALL_INC)
	     $(CPP) $(CPPFLAGS) -c $(OPPIC_PATH)/$(SRC)/core/oppic_util.cpp -o $(OBJ)/oppic_util.o $(ALL_INC)
	     $(CPP) $(CPPFLAGS) -c $(OPPIC_PATH)/$(SRC)/core/trace.cpp -o $(OBJ)/trace.o $(ALL_INC)
		 $(CPP) $(CPPFLAGS) -c $(OPPIC_PATH)/$(SRC)/core/opp_params.cpp -o $(OBJ)/opp_params.o $(ALL_INC)

oppic_seq: mklib oppic_core
	     $(CPP) $(CPPFLAGS) -c $(OPPIC_PATH)/$(SRC)/seq/oppic_seq.cpp -o $(OBJ)/oppic_seq.o $(ALL_INC)
	     ar -r $(LIB)/liboppic_seq.a $(OBJ)/oppic_lib_core.o $(OBJ)/oppic_util.o $(OBJ)/trace.o $(OBJ)/opp_params.o $(OBJ)/oppic_seq.o

oppic_omp: mklib oppic_core
	     $(CPP) $(CPPFLAGS) -c $(OPPIC_PATH)/$(SRC)/omp/oppic_omp.cpp -o $(OBJ)/oppic_omp.o $(ALL_INC)
	     ar -r $(LIB)/liboppic_omp.a $(OBJ)/oppic_lib_core.o $(OBJ)/oppic_util.o $(OBJ)/trace.o $(OBJ)/opp_params.o $(OBJ)/oppic_omp.o

oppic_cuda: mklib oppic_core
	     $(CPP) $(CPPFLAGS) -c $(OPPIC_PATH)/$(SRC)/cuda/oppic_cuda.cpp -o $(OBJ)/oppic_cuda.o $(ALL_INC)
		 nvcc $(NVCCFLAGS) -c $(OPPIC_PATH)/$(SRC)/cuda/particle_sorter.cu -o $(OBJ)/particle_sorter.o $(ALL_INC)
	     ar -r $(LIB)/liboppic_cuda.a $(OBJ)/oppic_lib_core.o $(OBJ)/oppic_util.o $(OBJ)/trace.o $(OBJ)/opp_params.o $(OBJ)/oppic_cuda.o $(OBJ)/particle_sorter.o

clean: 
	     rm -f *.o *.d *.a
	     rm -f $(OBJ)/*
	     rm -f $(LIB)/*

