
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
# scorep-mpiicpc  -qopt-report=2 -qopt-report-phase=vec -xHOST     --generate-line-info  -fsanitize=address
CPP        = mpicxx
MPICPP	   = mpicxx
NVCC       = nvcc
HIPCC      = hipcc

CPPFLAGS   = -std=c++14 -Wall -fPIC
NVCCFLAGS  = -m64 -Xptxas -dlcm=ca -Xptxas=-v $(NVCCFLAGS_ADD)
HIPCCFLAGS = -std=c++17 -isystem -D__HIP_PLATFORM_AMD__ $(NVCCFLAGS_ADD)

DEBUG ?= 0
ifeq ($(DEBUG), 1)
	CPPFLAGS   += -O0 -g
	NVCCFLAGS  += -O0 -g -lineinfo
	HIPCCFLAGS += -O0 -g -lineinfo
	DEBUG_LOG   = 1
else
	CPPFLAGS   += -O3 -g
	NVCCFLAGS  += -O3 -g
	HIPCCFLAGS += -O3 -g
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
	NVCCFLAGS += -DUSE_PETSC
	HIPCCFLAGS += -DUSE_PETSC
	PETSC_INC  = -I$(PETSC_INSTALL_PATH)/include
endif

CUDA ?= 0
ifeq ($(CUDA), 1)
	CUDA_INC   = -I$(CUDA_INSTALL_PATH)/include
endif

HIP ?= 1
ifeq ($(HIP), 1)
	HIP_INC   = -I$(ROCM_INSTALL_DIR)/include -I$(ROCM_THRUST_DIR)/include -I$(ROCM_PRIM_DIR)/include
endif

OPP_INC    = -I$(OPPIC_PATH)/include

PARMETIS ?= 0
ifeq ($(PARMETIS), 1)
	OPP_INC  += -DHAVE_PARMETIS
endif

.PHONY: clean mklib

ALL_TARGETS = clean seq omp cuda

all: $(ALL_TARGETS)

mklib:
	@mkdir -p $(LIB) $(OBJ)

seq: mklib
	     $(CPP) $(CPPFLAGS) -c $(SRC)/core/oppic_lib_core.cpp -o $(OBJ)/oppic_lib_core.o $(OPP_INC) $(PETSC_INC)
	     $(CPP) $(CPPFLAGS) -c $(SRC)/core/oppic_util.cpp -o $(OBJ)/oppic_util.o $(OPP_INC) $(PETSC_INC)
		 $(CPP) $(CPPFLAGS) -c $(SRC)/core/opp_params.cpp -o $(OBJ)/opp_params.o $(OPP_INC) $(PETSC_INC)
	     $(CPP) $(CPPFLAGS) -c $(SRC)/seq/oppic_seq.cpp -o $(OBJ)/oppic_seq.o $(OPP_INC) $(PETSC_INC)
	     ar -r $(LIB)/liboppic_seq.a $(OBJ)/oppic_lib_core.o $(OBJ)/oppic_util.o $(OBJ)/opp_params.o $(OBJ)/oppic_seq.o

omp: mklib
	     $(CPP) $(CPPFLAGS) -fopenmp -c $(SRC)/core/oppic_lib_core.cpp -o $(OBJ)/oppic_lib_core.o $(OPP_INC) $(PETSC_INC)
	     $(CPP) $(CPPFLAGS) -fopenmp -c $(SRC)/core/oppic_util.cpp -o $(OBJ)/oppic_util.o $(OPP_INC) $(PETSC_INC)
		 $(CPP) $(CPPFLAGS) -fopenmp -c $(SRC)/core/opp_params.cpp -o $(OBJ)/opp_params.o $(OPP_INC) $(PETSC_INC)
	     $(CPP) $(CPPFLAGS) -fopenmp -c $(SRC)/omp/oppic_omp.cpp -o $(OBJ)/oppic_omp.o $(OPP_INC) $(PETSC_INC)
	     ar -r $(LIB)/liboppic_omp.a $(OBJ)/oppic_lib_core.o $(OBJ)/oppic_util.o $(OBJ)/opp_params.o $(OBJ)/oppic_omp.o

cuda: mklib
	     $(CPP) $(CPPFLAGS) -DUSE_THRUST -c $(SRC)/core/oppic_lib_core.cpp -o $(OBJ)/oppic_lib_core.o $(OPP_INC) $(PETSC_INC) $(CUDA_INC)
	     $(CPP) $(CPPFLAGS) -DUSE_THRUST -c $(SRC)/core/oppic_util.cpp -o $(OBJ)/oppic_util.o $(OPP_INC) $(PETSC_INC) $(CUDA_INC)
		 $(CPP) $(CPPFLAGS) -DUSE_THRUST -c $(SRC)/core/opp_params.cpp -o $(OBJ)/opp_params.o $(OPP_INC) $(PETSC_INC) $(CUDA_INC)
	     $(NVCC) $(NVCCFLAGS) -DUSE_THRUST -c $(SRC)/cuda/opp_cuda.cu -o $(OBJ)/opp_cuda.o $(OPP_INC) $(PETSC_INC) $(CUDA_INC)
		 $(NVCC) $(NVCCFLAGS) -DUSE_THRUST -c $(SRC)/cuda/opp_particle_sorter.cu -o $(OBJ)/particle_sorter_cuda.o $(OPP_INC) $(PETSC_INC) $(CUDA_INC)
		 $(NVCC) $(NVCCFLAGS) -DUSE_THRUST -c $(SRC)/cuda/opp_cuda_halo.cu -o $(OBJ)/opp_cuda_halo.o $(OPP_INC) $(PETSC_INC) $(CUDA_INC)
	     ar -r $(LIB)/libopp_cuda.a $(OBJ)/oppic_lib_core.o $(OBJ)/oppic_util.o $(OBJ)/opp_params.o $(OBJ)/opp_cuda.o $(OBJ)/particle_sorter_cuda.o $(OBJ)/opp_cuda_halo.o

hip: mklib
	     $(HIPCC) $(HIPCCFLAGS) -DUSE_THRUST -c $(SRC)/core/oppic_lib_core.cpp -o $(OBJ)/oppic_lib_core.o $(OPP_INC) $(PETSC_INC) $(HIP_INC)
	     $(HIPCC) $(HIPCCFLAGS) -DUSE_THRUST -c $(SRC)/core/oppic_util.cpp -o $(OBJ)/oppic_util.o $(OPP_INC) $(PETSC_INC) $(HIP_INC)
		 $(HIPCC) $(HIPCCFLAGS) -DUSE_THRUST -c $(SRC)/core/opp_params.cpp -o $(OBJ)/opp_params.o $(OPP_INC) $(PETSC_INC) $(HIP_INC)
	     $(HIPCC) $(HIPCCFLAGS) -DUSE_THRUST -c $(SRC)/hip/opp_hip.cpp -o $(OBJ)/opp_hip.o $(OPP_INC) $(PETSC_INC) $(HIP_INC)
		 $(HIPCC) $(HIPCCFLAGS) -DUSE_THRUST -c $(SRC)/hip/opp_particle_sorter.cpp -o $(OBJ)/particle_sorter_hip.o $(OPP_INC) $(PETSC_INC) $(HIP_INC)
		 $(HIPCC) $(HIPCCFLAGS) -DUSE_THRUST -c $(SRC)/hip/opp_hip_halo.cpp -o $(OBJ)/opp_hip_halo.o $(OPP_INC) $(PETSC_INC) $(HIP_INC)
	     ar -r $(LIB)/libopp_hip.a $(OBJ)/oppic_lib_core.o $(OBJ)/oppic_util.o $(OBJ)/opp_params.o $(OBJ)/opp_hip.o $(OBJ)/particle_sorter_hip.o $(OBJ)/opp_hip_halo.o

mpi: mklib
	     $(CPP) $(CPPFLAGS) -DUSE_MPI -c $(SRC)/core/oppic_lib_core.cpp -o $(OBJ)/oppic_lib_core.o $(OPP_INC) $(PETSC_INC)
	     $(CPP) $(CPPFLAGS) -DUSE_MPI -c $(SRC)/core/oppic_util.cpp -o $(OBJ)/oppic_util.o $(OPP_INC) $(PETSC_INC)
		 $(CPP) $(CPPFLAGS) -DUSE_MPI -c $(SRC)/core/opp_params.cpp -o $(OBJ)/opp_params.o $(OPP_INC) $(PETSC_INC)
		 $(CPP) $(CPPFLAGS) -DUSE_MPI -c $(SRC)/mpi/opp_mpi_halo.cpp -o $(OBJ)/opp_mpi_halo.o $(OPP_INC) $(PETSC_INC)
		 $(CPP) $(CPPFLAGS) -DUSE_MPI -c $(SRC)/mpi/opp_mpi_halo_core.cpp -o $(OBJ)/opp_mpi_halo_core.o $(OPP_INC) $(PETSC_INC)
		 $(CPP) $(CPPFLAGS) -DUSE_MPI -c $(SRC)/mpi/opp_mpi_partition.cpp -o $(OBJ)/opp_mpi_partition.o $(OPP_INC) $(PETSC_INC)
		 $(CPP) $(CPPFLAGS) -DUSE_MPI -c $(SRC)/mpi/opp_mpi.cpp -o $(OBJ)/opp_mpi.o $(OPP_INC) $(PETSC_INC)
		 $(CPP) $(CPPFLAGS) -DUSE_MPI -c $(SRC)/mpi/opp_mpi_core.cpp -o $(OBJ)/opp_mpi_core.o $(OPP_INC) $(PETSC_INC)
		 $(CPP) $(CPPFLAGS) -DUSE_MPI -c $(SRC)/mpi/opp_mpi_utils.cpp -o $(OBJ)/opp_mpi_utils.o $(OPP_INC) $(PETSC_INC)
		 $(CPP) $(CPPFLAGS) -DUSE_MPI -c $(SRC)/mpi/opp_mpi_particle_comm.cpp -o $(OBJ)/opp_mpi_particle_comm.o $(OPP_INC) $(PETSC_INC)
		 $(CPP) $(CPPFLAGS) -DUSE_MPI -c $(SRC)/mpi/opp_mpi_double_ind_reducs.cpp -o $(OBJ)/opp_mpi_double_ind_reducs.o $(OPP_INC) $(PETSC_INC)
	     ar -r $(LIB)/liboppic_mpi.a $(OBJ)/oppic_lib_core.o $(OBJ)/oppic_util.o $(OBJ)/opp_params.o $(OBJ)/opp_mpi.o $(OBJ)/opp_mpi_core.o $(OBJ)/opp_mpi_utils.o $(OBJ)/opp_mpi_halo.o $(OBJ)/opp_mpi_halo_core.o $(OBJ)/opp_mpi_partition.o $(OBJ)/opp_mpi_particle_comm.o $(OBJ)/opp_mpi_double_ind_reducs.o

cuda_mpi: mklib
	     $(CPP) $(CPPFLAGS) -DUSE_MPI -DUSE_THRUST -c $(SRC)/core/oppic_lib_core.cpp -o $(OBJ)/oppic_lib_core.o $(OPP_INC) $(PETSC_INC) $(CUDA_INC)
	     $(CPP) $(CPPFLAGS) -DUSE_MPI -DUSE_THRUST -c $(SRC)/core/oppic_util.cpp -o $(OBJ)/oppic_util.o $(OPP_INC) $(PETSC_INC) $(CUDA_INC)
		 $(CPP) $(CPPFLAGS) -DUSE_MPI -DUSE_THRUST -c $(SRC)/core/opp_params.cpp -o $(OBJ)/opp_params.o $(OPP_INC) $(PETSC_INC) $(CUDA_INC)
		 $(CPP) $(CPPFLAGS) -DUSE_MPI -DUSE_THRUST -c $(SRC)/mpi/opp_mpi_halo_core.cpp -o $(OBJ)/opp_mpi_halo_core.o $(OPP_INC) $(PETSC_INC)
		 $(CPP) $(CPPFLAGS) -DUSE_MPI -DUSE_THRUST -c $(SRC)/mpi/opp_mpi_partition.cpp -o $(OBJ)/opp_mpi_partition.o $(OPP_INC) $(PETSC_INC)
		 $(CPP) $(CPPFLAGS) -DUSE_MPI -DUSE_THRUST -c $(SRC)/mpi/opp_mpi_core.cpp -o $(OBJ)/opp_mpi_core.o $(OPP_INC) $(PETSC_INC)
		 $(CPP) $(CPPFLAGS) -DUSE_MPI -DUSE_THRUST -c $(SRC)/mpi/opp_mpi_particle_comm.cpp -o $(OBJ)/opp_mpi_particle_comm.o $(OPP_INC) $(PETSC_INC)
		 $(CPP) $(CPPFLAGS) -DUSE_MPI -DUSE_THRUST -c $(SRC)/mpi/opp_mpi_double_ind_reducs.cpp -o $(OBJ)/opp_mpi_double_ind_reducs.o $(OPP_INC) $(PETSC_INC)
		 $(CPP) $(CPPFLAGS) -DUSE_MPI -DUSE_THRUST -c $(SRC)/mpi/opp_mpi_utils.cpp -o $(OBJ)/opp_mpi_utils.o $(OPP_INC) $(PETSC_INC)
	     $(NVCC) $(NVCCFLAGS) -DUSE_MPI -DUSE_THRUST -c $(SRC)/cuda/opp_cuda.cu -o $(OBJ)/opp_cuda.o $(OPP_INC) $(PETSC_INC) $(CUDA_INC)
		 $(NVCC) $(NVCCFLAGS) -DUSE_MPI -DUSE_THRUST -c $(SRC)/cuda/opp_cuda_halo.cu -o $(OBJ)/opp_cuda_halo.o $(OPP_INC) $(PETSC_INC) $(CUDA_INC)
		 $(NVCC) $(NVCCFLAGS) -DUSE_MPI -DUSE_THRUST -c $(SRC)/cuda/opp_particle_sorter.cu -o $(OBJ)/particle_sorter_cuda.o $(OPP_INC) $(PETSC_INC) $(CUDA_INC)
	     ar -r $(LIB)/libopp_cuda_mpi.a $(OBJ)/oppic_lib_core.o $(OBJ)/oppic_util.o $(OBJ)/opp_params.o $(OBJ)/opp_mpi_halo_core.o $(OBJ)/opp_mpi_partition.o $(OBJ)/opp_mpi_core.o $(OBJ)/opp_mpi_particle_comm.o $(OBJ)/opp_mpi_double_ind_reducs.o $(OBJ)/opp_mpi_utils.o $(OBJ)/opp_cuda.o $(OBJ)/opp_cuda_halo.o $(OBJ)/particle_sorter_cuda.o				 
		 
clean: 
	     rm -f *.o *.d *.a
	     rm -f $(OBJ)/*
	     rm -f $(LIB)/*

