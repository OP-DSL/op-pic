
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
NVCC       = nvcc
HIPCC      = hipcc

CPPFLAGS   = -std=c++11 -Wall -fopenmp -fPIC
NVCCFLAGS  = -m64 -Xptxas -dlcm=ca -Xptxas=-v $(NVCCFLAGS_ADD)
HIPCCFLAGS = -std=c++11 -isystem -D__HIP_PLATFORM_NVIDIA__ $(NVCCFLAGS_ADD)

DEBUG ?= 0
ifeq ($(DEBUG), 1)
	CPPFLAGS   += -O0 -g -debug all
	NVCCFLAGS  += -O0 -g -lineinfo
	HIPCCFLAGS += -O0 -g -lineinfo
	DEBUG_LOG   = 1
else
	CPPFLAGS   += -O3 
	NVCCFLAGS  += -O3
	HIPCCFLAGS += -O3
endif

IEEE ?= 0
ifeq ($(IEEE), 1)
	CPPFLAGS   += -fp-model strict -fp-model source -prec-div -prec-sqrt
	NVCCFLAGS  += -fmad=false -ftz=false -prec-div=true -prec-sqrt=true
	HIPCCFLAGS  += -fmad=false -ftz=false -prec-div=true -prec-sqrt=true
endif

DEBUG_LOG ?= 0
ifeq ($(DEBUG_LOG), 1)
	CPPFLAGS   += -DDEBUG_LOG
	NVCCFLAGS  += -DDEBUG_LOG
	HIPCCFLAGS += -DDEBUG_LOG
endif

PETSC ?= 1
ifeq ($(PETSC), 1)
	CPPFLAGS   += -DUSE_PETSC
	# HIPCCFLAGS += -DUSE_PETSC
	PETSC_INC  = -I$(PETSC_INSTALL_PATH)/include
endif

CUDA ?= 1
ifeq ($(CUDA), 1)
	CUDA_INC   = -I$(CUDA_INSTALL_PATH)/include
endif

HIP ?= 1
ifeq ($(HIP), 1)
	HIP_INC   = -I$(ROCM_INSTALL_DIR)/include -I$(ROCM_THRUST_DIR)/include -I$(ROCM_PRIM_DIR)/include
endif

OPP_INC    = -I$(OPPIC_PATH)/include

.PHONY: clean mklib

ALL_TARGETS = clean seq omp cuda

all: $(ALL_TARGETS)

mklib:
	@mkdir -p $(LIB) $(OBJ)

seq: mklib
	     $(CPP) $(CPPFLAGS) -c $(SRC)/core/oppic_lib_core.cpp -o $(OBJ)/oppic_lib_core.o $(OPP_INC) $(PETSC_INC)
	     $(CPP) $(CPPFLAGS) -c $(SRC)/core/oppic_util.cpp -o $(OBJ)/oppic_util.o $(OPP_INC) $(PETSC_INC)
	     $(CPP) $(CPPFLAGS) -c $(SRC)/core/trace.cpp -o $(OBJ)/trace.o $(OPP_INC) $(PETSC_INC)
		 $(CPP) $(CPPFLAGS) -c $(SRC)/core/opp_params.cpp -o $(OBJ)/opp_params.o $(OPP_INC) $(PETSC_INC)
	     $(CPP) $(CPPFLAGS) -c $(SRC)/seq/oppic_seq.cpp -o $(OBJ)/oppic_seq.o $(OPP_INC) $(PETSC_INC)
	     ar -r $(LIB)/liboppic_seq.a $(OBJ)/oppic_lib_core.o $(OBJ)/oppic_util.o $(OBJ)/trace.o $(OBJ)/opp_params.o $(OBJ)/oppic_seq.o

omp: mklib
	     $(CPP) $(CPPFLAGS) -c $(SRC)/core/oppic_lib_core.cpp -o $(OBJ)/oppic_lib_core.o $(OPP_INC) $(PETSC_INC)
	     $(CPP) $(CPPFLAGS) -c $(SRC)/core/oppic_util.cpp -o $(OBJ)/oppic_util.o $(OPP_INC) $(PETSC_INC)
	     $(CPP) $(CPPFLAGS) -c $(SRC)/core/trace.cpp -o $(OBJ)/trace.o $(OPP_INC) $(PETSC_INC)
		 $(CPP) $(CPPFLAGS) -c $(SRC)/core/opp_params.cpp -o $(OBJ)/opp_params.o $(OPP_INC) $(PETSC_INC)
	     $(CPP) $(CPPFLAGS) -c $(SRC)/omp/oppic_omp.cpp -o $(OBJ)/oppic_omp.o $(OPP_INC) $(PETSC_INC)
	     ar -r $(LIB)/liboppic_omp.a $(OBJ)/oppic_lib_core.o $(OBJ)/oppic_util.o $(OBJ)/trace.o $(OBJ)/opp_params.o $(OBJ)/oppic_omp.o

cuda: mklib
	     $(CPP) $(CPPFLAGS) -DUSE_THRUST -c $(SRC)/core/oppic_lib_core.cpp -o $(OBJ)/oppic_lib_core.o $(OPP_INC) $(PETSC_INC) $(CUDA_INC)
	     $(CPP) $(CPPFLAGS) -DUSE_THRUST -c $(SRC)/core/oppic_util.cpp -o $(OBJ)/oppic_util.o $(OPP_INC) $(PETSC_INC) $(CUDA_INC)
	     $(CPP) $(CPPFLAGS) -DUSE_THRUST -c $(SRC)/core/trace.cpp -o $(OBJ)/trace.o $(OPP_INC) $(PETSC_INC) $(CUDA_INC)
		 $(CPP) $(CPPFLAGS) -DUSE_THRUST -c $(SRC)/core/opp_params.cpp -o $(OBJ)/opp_params.o $(OPP_INC) $(PETSC_INC) $(CUDA_INC)
	     $(NVCC) $(NVCCFLAGS) -DUSE_THRUST -c $(SRC)/cuda/oppic_cuda.cu -o $(OBJ)/oppic_cuda.o $(OPP_INC) $(PETSC_INC) $(CUDA_INC)
		 $(NVCC) $(NVCCFLAGS) -DUSE_THRUST -c $(SRC)/cuda/particle_sorter.cu -o $(OBJ)/particle_sorter_cuda.o $(OPP_INC) $(PETSC_INC) $(CUDA_INC)
	     ar -r $(LIB)/liboppic_cuda.a $(OBJ)/oppic_lib_core.o $(OBJ)/oppic_util.o $(OBJ)/trace.o $(OBJ)/opp_params.o $(OBJ)/oppic_cuda.o $(OBJ)/particle_sorter_cuda.o

hip: mklib
	     $(HIPCC) $(HIPCCFLAGS) -DUSE_THRUST -c $(SRC)/hip/oppic_hip.cpp -o $(OBJ)/oppic_hip.o $(CUDA_INC) $(HIP_INC) $(OPP_INC) $(PETSC_INC)
		 $(HIPCC) $(HIPCCFLAGS) -DUSE_THRUST -c $(SRC)/hip/particle_sorter.cpp -o $(OBJ)/particle_sorter_hip.o $(HIP_INC) $(OPP_INC) $(PETSC_INC) $(CUDA_INC)
	     $(HIPCC) $(HIPCCFLAGS) -DUSE_THRUST -c $(SRC)/core/oppic_lib_core.cpp -o $(OBJ)/oppic_lib_core.o $(OPP_INC) $(PETSC_INC) $(HIP_INC) $(CUDA_INC)
	     $(HIPCC) $(HIPCCFLAGS) -DUSE_THRUST -c $(SRC)/core/oppic_util.cpp -o $(OBJ)/oppic_util.o $(OPP_INC) $(PETSC_INC) $(HIP_INC) $(CUDA_INC)
	     $(HIPCC) $(HIPCCFLAGS) -DUSE_THRUST -c $(SRC)/core/trace.cpp -o $(OBJ)/trace.o $(OPP_INC) $(PETSC_INC) $(HIP_INC) $(CUDA_INC)
		 $(HIPCC) $(HIPCCFLAGS) -DUSE_THRUST -c $(SRC)/core/opp_params.cpp -o $(OBJ)/opp_params.o $(OPP_INC) $(PETSC_INC) $(HIP_INC) $(CUDA_INC)
	     ar -r $(LIB)/liboppic_hip.a $(OBJ)/oppic_lib_core.o $(OBJ)/oppic_util.o $(OBJ)/trace.o $(OBJ)/opp_params.o $(OBJ)/oppic_hip.o $(OBJ)/particle_sorter_hip.o

mpi: mklib
	     $(CPP) $(CPPFLAGS) -c $(SRC)/core/oppic_lib_core.cpp -o $(OBJ)/oppic_lib_core.o $(OPP_INC) $(PETSC_INC)
	     $(CPP) $(CPPFLAGS) -c $(SRC)/core/oppic_util.cpp -o $(OBJ)/oppic_util.o $(OPP_INC) $(PETSC_INC)
	     $(CPP) $(CPPFLAGS) -c $(SRC)/core/trace.cpp -o $(OBJ)/trace.o $(OPP_INC) $(PETSC_INC)
		 $(CPP) $(CPPFLAGS) -c $(SRC)/core/opp_params.cpp -o $(OBJ)/opp_params.o $(OPP_INC) $(PETSC_INC)
		 $(CPP) $(CPPFLAGS) -c $(SRC)/mpi/opp_mpi_halo.cpp -o $(OBJ)/opp_mpi_halo.o $(OPP_INC) $(PETSC_INC)
		 $(CPP) $(CPPFLAGS) -c $(SRC)/mpi/opp_mpi_partition.cpp -o $(OBJ)/opp_mpi_partition.o $(OPP_INC) $(PETSC_INC)
		 $(CPP) $(CPPFLAGS) -c $(SRC)/mpi/opp_mpi.cpp -o $(OBJ)/opp_mpi.o $(OPP_INC) $(PETSC_INC)
		 $(CPP) $(CPPFLAGS) -c $(SRC)/mpi/opp_mpi_particle_comm.cpp -o $(OBJ)/opp_mpi_particle_comm.o $(OPP_INC) $(PETSC_INC)
	     ar -r $(LIB)/liboppic_mpi.a $(OBJ)/oppic_lib_core.o $(OBJ)/oppic_util.o $(OBJ)/trace.o $(OBJ)/opp_params.o $(OBJ)/opp_mpi.o $(OBJ)/opp_mpi_halo.o $(OBJ)/opp_mpi_partition.o $(OBJ)/opp_mpi_particle_comm.o
		 
		 
clean: 
	     rm -f *.o *.d *.a
	     rm -f $(OBJ)/*
	     rm -f $(LIB)/*

