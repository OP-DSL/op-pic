
# source /ext-home/zl/phd/OP-PIC/scripts/source_oneapi

SRC = .
INC = .
LIB = lib
OBJ = obj
BIN = bin
# scorep-mpiicpc  -qopt-report=2 -qopt-report-phase=vec -xHOST          --generate-line-info
CPP        = mpicxx
MPICPP	   = mpicxx
NVCC       = nvcc
HIPCC      = hipcc

CPPFLAGS   = -std=c++17 -Wall -fPIC
NVCCFLAGS  = -m64 -Xptxas -dlcm=ca -Xptxas=-v $(NVCCFLAGS_ADD)
HIPCCFLAGS = -std=c++17 -Wall -isystem -D__HIP_PLATFORM_AMD__

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

# TODO: Remove; Only for development
CPU_ASSIGN ?= 0
ifeq ($(CPU_ASSIGN), 1)
	NVCCFLAGS  += -DUSE_CPU_ASSIGN
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
	PETSC_INC  = -I$(PETSC_INSTALL_PATH)/include
	PETSC_LIB  = -L$(PETSC_INSTALL_PATH)/lib -lpetsc
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

PARMETIS ?= 0
ifeq ($(PARMETIS), 1)
	ALL_LIBS = -lparmetis -lmetis
	ALL_INC  += -DHAVE_PARMETIS
endif

ALL_INC    += -I$(OPPIC_PATH)/include $(PETSC_INC) -I/ext-home/zl/lib/inclide
ALL_LIBS   += -L$(OPPIC_PATH)/lib -L/ext-home/zl/lib/lib $(PETSC_LIB) 



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

t:
	     $(CPP) $(CPPFLAGS) -c fempic_test.cpp -o $(OBJ)/fempic_test.o $(ALL_INC)
	     $(CPP) $(CPPFLAGS) -c fempic_ori/meshes.cpp -o $(OBJ)/meshes.o $(ALL_INC)
	     $(CPP) $(CPPFLAGS) -o $(BIN)/test $(OBJ)/fempic_test.o $(OBJ)/meshes.o $(OBJ)/particles.o $(OBJ)/maths.o $(OBJ)/FESolver.o $(ALL_INC) $(ALL_LIBS) -loppic_seq


omp: common
		 $(CPP) $(CPPFLAGS) -fopenmp -c omp/ompkernels.cpp -o $(OBJ)/ompkernels.o $(ALL_INC)
	     $(CPP) $(CPPFLAGS) -fopenmp -o $(BIN)/omp $(OBJ)/fempic_op.o $(OBJ)/meshes.o $(OBJ)/particles.o $(OBJ)/maths.o $(OBJ)/ompkernels.o $(OBJ)/FESolver.o $(ALL_INC) $(ALL_LIBS) -loppic_omp

mpiX: common
		 $(CPP) $(CPPFLAGS) -c mpi/mpikernels.cpp -o $(OBJ)/mpikernels.o $(ALL_INC)
	     $(CPP) $(CPPFLAGS) -o $(BIN)/mpi $(OBJ)/fempic_op.o $(OBJ)/meshes.o $(OBJ)/particles.o $(OBJ)/maths.o $(OBJ)/mpikernels.o $(OBJ)/FESolver.o  -loppic_mpi $(ALL_INC) $(ALL_LIBS)

mpi_test: mklib
	     $(CPP) $(CPPFLAGS) -c fempic_test.cpp -o $(OBJ)/fempic_test.o $(ALL_INC)
		 $(CPP) $(CPPFLAGS) -c fempic_ori/meshes.cpp -o $(OBJ)/meshes.o $(ALL_INC)
	     $(CPP) $(CPPFLAGS) -c fempic_ori/particles.cpp -o $(OBJ)/particles.o $(ALL_INC)
		 $(CPP) $(CPPFLAGS) -c fempic_ori/maths.cpp -o $(OBJ)/maths.o $(ALL_INC)
		 $(CPP) $(CPPFLAGS) -c FESolver.cpp -o $(OBJ)/FESolver.o $(ALL_INC)		 
		 $(CPP) $(CPPFLAGS) -c mpi/mpikernels.cpp -o $(OBJ)/mpikernels.o $(ALL_INC)
	     $(CPP) $(CPPFLAGS) -o $(BIN)/mpi $(OBJ)/fempic_test.o $(OBJ)/meshes.o $(OBJ)/particles.o $(OBJ)/maths.o $(OBJ)/mpikernels.o $(OBJ)/FESolver.o  -loppic_mpi $(ALL_INC) $(ALL_LIBS)

mpi: mklib
	     $(CPP) $(CPPFLAGS) -DENABLE_MPI -c fempic_op.cpp -o $(OBJ)/fempic_op.o $(ALL_INC)
		 $(CPP) $(CPPFLAGS) -DENABLE_MPI -c fempic_ori/meshes.cpp -o $(OBJ)/meshes.o $(ALL_INC)
	     $(CPP) $(CPPFLAGS) -DENABLE_MPI -c fempic_ori/particles.cpp -o $(OBJ)/particles.o $(ALL_INC)
		 $(CPP) $(CPPFLAGS) -DENABLE_MPI -c fempic_ori/maths.cpp -o $(OBJ)/maths.o $(ALL_INC)
		 $(CPP) $(CPPFLAGS) -DENABLE_MPI -c FESolver.cpp -o $(OBJ)/FESolver.o $(ALL_INC)		 
		 $(CPP) $(CPPFLAGS) -DENABLE_MPI -c mpi/mpikernels.cpp -o $(OBJ)/mpikernels.o $(ALL_INC)
	     $(CPP) $(CPPFLAGS) -DENABLE_MPI -o $(BIN)/mpi $(OBJ)/fempic_op.o $(OBJ)/meshes.o $(OBJ)/particles.o $(OBJ)/maths.o $(OBJ)/mpikernels.o $(OBJ)/FESolver.o  -loppic_mpi $(ALL_INC) $(ALL_LIBS)

cuda: common
		 $(NVCC) $(NVCCFLAGS) -c cuda/cudakernels.cu -o $(OBJ)/cudakernels.o $(ALL_INC) $(CUDA_INC)
	     $(CPP) $(CPPFLAGS) -o $(BIN)/cuda $(OBJ)/fempic_op.o $(OBJ)/meshes.o $(OBJ)/particles.o $(OBJ)/maths.o $(OBJ)/cudakernels.o $(OBJ)/FESolver.o $(ALL_INC) $(CUDA_INC) $(ALL_LIBS) $(CUDA_LIB) -loppic_cuda

hip: mklib
		 $(HIPCC) $(HIPCCFLAGS) -c hip/hipkernels.cpp -o $(OBJ)/hipkernels.o $(ALL_INC) $(HIP_INC)
	     $(HIPCC) $(HIPCCFLAGS) -c fempic_op.cpp -o $(OBJ)/fempic_op.o $(ALL_INC)
	     $(HIPCC) $(HIPCCFLAGS) -c fempic_ori/meshes.cpp -o $(OBJ)/meshes.o $(ALL_INC)
	     $(HIPCC) $(HIPCCFLAGS) -c fempic_ori/particles.cpp -o $(OBJ)/particles.o $(ALL_INC)
		 $(HIPCC) $(HIPCCFLAGS) -c fempic_ori/maths.cpp -o $(OBJ)/maths.o $(ALL_INC)
		 $(HIPCC) $(HIPCCFLAGS) -c FESolver.cpp -o $(OBJ)/FESolver.o $(ALL_INC)
	     $(HIPCC) $(HIPCCFLAGS) -o $(BIN)/hip $(OBJ)/fempic_op.o $(OBJ)/meshes.o $(OBJ)/particles.o $(OBJ)/maths.o $(OBJ)/hipkernels.o $(OBJ)/FESolver.o $(ALL_INC) $(HIP_INC) $(ALL_LIBS) $(HIP_LIB) -loppic_hip

clean: 
	     rm -f *.o *.d *.a
	     rm -f $(OBJ)/*
	     rm -f $(BIN)/*
	     rm -f $(LIB)/*
		 
