# export OPPIC_PATH=/ext-home/zl/phd/neso_test/OP-PIC/lib_oppic

SRC = .
INC = .
LIB = lib
OBJ = obj
BIN = bin

CPP        = $(CC_COMPILER)
MPICPP	   = $(MPI_COMPILER)
NVCC       = nvcc
HIPCC      = $(HIP_COMPILER)
SYCLCC     = icpx

CPPFLAGS   = -std=c++14 -Wall -fPIC 
NVCCFLAGS  = -m64 -Xptxas -dlcm=ca -Xptxas=-v $(NVCCFLAGS_ADD)
HIPCCFLAGS = $(HIPCCFLAGS_ADD) -std=c++17 -Wall -isystem -D__HIP_PLATFORM_AMD__
SYCLFLAGS = -w -fsycl -std=c++17 -Wall -fPIC -march=native 

DEBUG ?= 0
ifeq ($(DEBUG), 1)
	CPPFLAGS   += -O0 -g
	NVCCFLAGS  += -O0 -g -lineinfo
	HIPCCFLAGS += -O0 -g
	DEBUG_LOG   = 1
else
	CPPFLAGS   += -O3 -g
	NVCCFLAGS  += -O3 -g
	HIPCCFLAGS += -O3 -g
endif

DEBUG_LOG ?= 0
ifeq ($(DEBUG_LOG), 1)
	CPPFLAGS   += -DDEBUG_LOG
	NVCCFLAGS  += -DDEBUG_LOG
	HIPCCFLAGS += -DDEBUG_LOG
endif

IEEE ?= 0
ifeq ($(IEEE), 1)
	CPPFLAGS   += -fp-model strict -fp-model source -prec-div -prec-sqrt
	NVCCFLAGS  += -fmad=false -ftz=false -prec-div=true -prec-sqrt=true
	HIPCCFLAGS  += -fmad=false -ftz=false -prec-div=true -prec-sqrt=true
endif

CUDA ?= 1
ifeq ($(CUDA), 1)
	CUDA_INC   = -I$(CUDA_INSTALL_PATH)/include
	CUDA_LIB   = -L$(CUDA_INSTALL_PATH)/lib64 -lcudart -lcuda
endif

HIP ?= 1
ifeq ($(HIP), 1)
	HIP_INC   = -I$(ROCM_INSTALL_DIR)/include -I$(ROCM_THRUST_DIR)/include -I$(ROCM_PRIM_DIR)/include
endif

ALL_INC    += -I$(OPP_PATH)/include
ALL_LIBS   += -L$(OPP_PATH)/lib

.PHONY: clean mklib

ALL_TARGETS = clean seq omp cuda mpi

all: $(ALL_TARGETS)

mklib:
	@mkdir -p $(OBJ) $(BIN)

common: mklib 

seq: common
		 $(CPP) $(CPPFLAGS) -c advec_opp.cpp -o $(OBJ)/advec_opp_seq.o $(ALL_INC)
		 $(CPP) $(CPPFLAGS) -c seq/opp_kernels.cpp -o $(OBJ)/opp_kernels_seq.o $(ALL_INC)
	     $(CPP) $(CPPFLAGS) -o $(BIN)/seq \
		 	$(OBJ)/advec_opp_seq.o \
			$(OBJ)/opp_kernels_seq.o \
			$(ALL_INC) $(ALL_LIBS) -lopp_seq

mpi: common
		 $(MPICPP) $(CPPFLAGS) -DUSE_MPI -c advec_opp.cpp -o $(OBJ)/advec_opp_mpi.o $(ALL_INC)
		 $(MPICPP) $(CPPFLAGS) -DUSE_MPI -c mpi/opp_kernels.cpp -o $(OBJ)/opp_kernels_mpi.o $(ALL_INC)
	     $(MPICPP) $(CPPFLAGS) -DUSE_MPI -o $(BIN)/mpi \
		 	$(OBJ)/advec_opp_mpi.o \
			$(OBJ)/opp_kernels_mpi.o \
			$(ALL_INC) $(ALL_LIBS) -lopp_mpi

omp: common
	     $(CPP) $(CPPFLAGS) -fopenmp -DUSE_OMP -c advec_opp.cpp -o $(OBJ)/advec_opp_omp.o $(ALL_INC)
	     $(CPP) $(CPPFLAGS) -fopenmp -DUSE_OMP -c omp/opp_kernels.cpp -o $(OBJ)/opp_kernels_omp.o $(ALL_INC)
	     $(CPP) $(CPPFLAGS) -fopenmp -o $(BIN)/omp \
		 	$(OBJ)/advec_opp_omp.o \
			$(OBJ)/opp_kernels_omp.o \
			$(ALL_INC) $(ALL_LIBS) -lopp_omp

omp_mpi: common
	     $(CPP) $(CPPFLAGS) -fopenmp -DUSE_MPI -DUSE_OMP -c advec_opp.cpp -o $(OBJ)/advec_opp_omp_mpi.o $(ALL_INC)
	     $(CPP) $(CPPFLAGS) -fopenmp -DUSE_MPI -DUSE_OMP -c omp/opp_kernels.cpp -o $(OBJ)/opp_kernels_omp_mpi.o $(ALL_INC)
	     $(CPP) $(CPPFLAGS) -fopenmp -o $(BIN)/omp_mpi \
		 	$(OBJ)/advec_opp_omp_mpi.o \
			$(OBJ)/opp_kernels_omp_mpi.o \
			$(ALL_INC) $(ALL_LIBS) -lopp_omp_mpi

cuda: common
		 $(CPP) $(CPPFLAGS) -c advec_opp.cpp -o $(OBJ)/advec_opp_cuda.o $(ALL_INC)
		 $(NVCC) $(NVCCFLAGS) -c cuda/opp_kernels.cu -o $(OBJ)/opp_kernels_cuda.o $(ALL_INC) $(CUDA_INC)
	     $(CPP) $(CPPFLAGS) -o $(BIN)/cuda \
		 	$(OBJ)/advec_opp_cuda.o \
			$(OBJ)/opp_kernels_cuda.o \
			$(ALL_INC) $(CUDA_INC) $(ALL_LIBS) $(CUDA_LIB) -lopp_cuda

cuda_mpi: common
		 $(MPICPP) $(CPPFLAGS) -DUSE_MPI -c advec_opp.cpp -o $(OBJ)/advec_opp_cuda_mpi.o $(ALL_INC) $(CUDA_INC)
		 $(NVCC) $(NVCCFLAGS) -DUSE_MPI -c cuda/opp_kernels.cu -o $(OBJ)/opp_kernels_cuda_mpi.o $(ALL_INC) $(CUDA_INC)
	     $(MPICPP) $(CPPFLAGS) -DUSE_MPI -o $(BIN)/cuda_mpi \
		 	$(OBJ)/advec_opp_cuda_mpi.o \
			$(OBJ)/opp_kernels_cuda_mpi.o \
			$(ALL_INC) $(CUDA_INC) $(ALL_LIBS) $(CUDA_LIB) -lopp_cuda_mpi

sycl: common
		 $(SYCLCC) $(SYCLFLAGS) -c advec_opp.cpp -o $(OBJ)/advec_opp_cuda.o $(ALL_INC)
		 $(SYCLCC) $(SYCLFLAGS) -c sycl/opp_kernels.dp.cpp -o $(OBJ)/opp_kernels_cuda.o $(ALL_INC) $(CUDA_INC)
	     $(SYCLCC) $(SYCLFLAGS) -o $(BIN)/cuda \
		 	$(OBJ)/advec_opp_cuda.o \
			$(OBJ)/opp_kernels_cuda.o \
			$(ALL_INC) $(CUDA_INC) $(ALL_LIBS) $(CUDA_LIB) -lopp_sycl

sycl_mpi: common
		 $(SYCLCC) $(SYCLFLAGS) -DUSE_MPI -c advec_opp.cpp -o $(OBJ)/advec_opp_cuda_mpi.o $(ALL_INC) $(CUDA_INC)
		 $(SYCLCC) $(SYCLFLAGS) -DUSE_MPI -c sycl/opp_kernels.dp.cpp -o $(OBJ)/opp_kernels_cuda_mpi.o $(ALL_INC) $(CUDA_INC)
	     $(SYCLCC) $(SYCLFLAGS) -DUSE_MPI -o $(BIN)/cuda_mpi \
		 	$(OBJ)/advec_opp_cuda_mpi.o \
			$(OBJ)/opp_kernels_cuda_mpi.o \
			$(ALL_INC) $(CUDA_INC) $(ALL_LIBS) $(CUDA_LIB) -lopp_sycl_mpi -lmpi

hip: mklib
		 $(HIPCC) $(HIPCCFLAGS) -c advec_opp.cpp -o $(OBJ)/advec_opp_hip.o $(ALL_INC) $(HIP_INC)
		 $(HIPCC) $(HIPCCFLAGS) -c hip/opp_kernels.cpp -o $(OBJ)/opp_kernels_hip.o $(ALL_INC) $(HIP_INC)
	     $(HIPCC) $(HIPCCFLAGS) -o $(BIN)/hip \
			$(OBJ)/advec_opp_hip.o \
			$(OBJ)/opp_kernels_hip.o \
			$(ALL_INC) $(HIP_INC) $(ALL_LIBS) -lopp_hip

hip_mpi: common
		 $(HIPCC) $(HIPCCFLAGS) -DUSE_MPI -c advec_opp.cpp -o $(OBJ)/advec_opp_hip_mpi.o $(ALL_INC) $(HIP_INC)
		 $(HIPCC) $(HIPCCFLAGS) -DUSE_MPI -c hip/opp_kernels.cpp -o $(OBJ)/opp_kernels_hip_mpi.o $(ALL_INC) $(HIP_INC)
	     $(HIPCC) $(HIPCCFLAGS) -DUSE_MPI -o $(BIN)/hip_mpi \
		 	$(OBJ)/advec_opp_hip_mpi.o \
			$(OBJ)/opp_kernels_hip_mpi.o \
			$(ALL_INC) $(HIP_INC) $(ALL_LIBS) -lopp_hip_mpi -lmpi

clean: 
	     rm -f *.o *.d *.a
	     rm -f $(OBJ)/*
	     rm -f $(BIN)/*
	     rm -f $(LIB)/*
		 
