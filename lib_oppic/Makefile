
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
endif

CUDA ?= 1
ifeq ($(CUDA), 1)
	CUDA_INC   = -I$(CUDA_INSTALL_PATH)/include
endif

OPP_INC    = -I$(OPPIC_PATH)/include

.PHONY: clean mklib

ALL_TARGETS = clean oppic_seq oppic_omp oppic_cuda

all: $(ALL_TARGETS)

mklib:
	@mkdir -p $(LIB) $(OBJ)

oppic_seq: mklib
	     $(CPP) $(CPPFLAGS) -c $(SRC)/core/oppic_lib_core.cpp -o $(OBJ)/oppic_lib_core.o $(OPP_INC) $(PETSC_INC)
	     $(CPP) $(CPPFLAGS) -c $(SRC)/core/oppic_util.cpp -o $(OBJ)/oppic_util.o $(OPP_INC) $(PETSC_INC)
	     $(CPP) $(CPPFLAGS) -c $(SRC)/core/trace.cpp -o $(OBJ)/trace.o $(OPP_INC) $(PETSC_INC)
		 $(CPP) $(CPPFLAGS) -c $(SRC)/core/opp_params.cpp -o $(OBJ)/opp_params.o $(OPP_INC) $(PETSC_INC)
	     $(CPP) $(CPPFLAGS) -c $(SRC)/seq/oppic_seq.cpp -o $(OBJ)/oppic_seq.o $(OPP_INC) $(PETSC_INC)
	     ar -r $(LIB)/liboppic_seq.a $(OBJ)/oppic_lib_core.o $(OBJ)/oppic_util.o $(OBJ)/trace.o $(OBJ)/opp_params.o $(OBJ)/oppic_seq.o

oppic_omp: mklib
	     $(CPP) $(CPPFLAGS) -c $(SRC)/core/oppic_lib_core.cpp -o $(OBJ)/oppic_lib_core.o $(OPP_INC) $(PETSC_INC)
	     $(CPP) $(CPPFLAGS) -c $(SRC)/core/oppic_util.cpp -o $(OBJ)/oppic_util.o $(OPP_INC) $(PETSC_INC)
	     $(CPP) $(CPPFLAGS) -c $(SRC)/core/trace.cpp -o $(OBJ)/trace.o $(OPP_INC) $(PETSC_INC)
		 $(CPP) $(CPPFLAGS) -c $(SRC)/core/opp_params.cpp -o $(OBJ)/opp_params.o $(OPP_INC) $(PETSC_INC)
	     $(CPP) $(CPPFLAGS) -c $(SRC)/omp/oppic_omp.cpp -o $(OBJ)/oppic_omp.o $(OPP_INC) $(PETSC_INC)
	     ar -r $(LIB)/liboppic_omp.a $(OBJ)/oppic_lib_core.o $(OBJ)/oppic_util.o $(OBJ)/trace.o $(OBJ)/opp_params.o $(OBJ)/oppic_omp.o

oppic_cuda: mklib
	     $(CPP) $(CPPFLAGS) -DUSE_THRUST -c $(SRC)/core/oppic_lib_core.cpp -o $(OBJ)/oppic_lib_core.o $(OPP_INC) $(PETSC_INC) $(CUDA_INC)
	     $(CPP) $(CPPFLAGS) -DUSE_THRUST -c $(SRC)/core/oppic_util.cpp -o $(OBJ)/oppic_util.o $(OPP_INC) $(PETSC_INC) $(CUDA_INC)
	     $(CPP) $(CPPFLAGS) -DUSE_THRUST -c $(SRC)/core/trace.cpp -o $(OBJ)/trace.o $(OPP_INC) $(PETSC_INC) $(CUDA_INC)
		 $(CPP) $(CPPFLAGS) -DUSE_THRUST -c $(SRC)/core/opp_params.cpp -o $(OBJ)/opp_params.o $(OPP_INC) $(PETSC_INC) $(CUDA_INC)
	     nvcc $(NVCCFLAGS) -DUSE_THRUST -c $(SRC)/cuda/oppic_cuda.cu -o $(OBJ)/oppic_cuda.o $(OPP_INC) $(PETSC_INC) $(CUDA_INC)
		 nvcc $(NVCCFLAGS) -DUSE_THRUST -c $(SRC)/cuda/particle_sorter.cu -o $(OBJ)/particle_sorter.o $(OPP_INC) $(PETSC_INC) $(CUDA_INC)
	     ar -r $(LIB)/liboppic_cuda.a $(OBJ)/oppic_lib_core.o $(OBJ)/oppic_util.o $(OBJ)/trace.o $(OBJ)/opp_params.o $(OBJ)/oppic_cuda.o $(OBJ)/particle_sorter.o

clean: 
	     rm -f *.o *.d *.a
	     rm -f $(OBJ)/*
	     rm -f $(LIB)/*

