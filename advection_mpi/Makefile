
SRC = .
INC = .
LIB = lib
OBJ = obj
BIN = bin

CPP        = mpicxx
MPICPP	   = mpicxx
NVCC       = nvcc

CPPFLAGS   = -std=c++14 -Wall -fPIC 
NVCCFLAGS  = -m64 -Xptxas -dlcm=ca -Xptxas=-v $(NVCCFLAGS_ADD)

DEBUG ?= 0
ifeq ($(DEBUG), 1)
	CPPFLAGS   += -O0 -g
	NVCCFLAGS  += -O0 -g -lineinfo
	DEBUG_LOG   = 1
else
	CPPFLAGS   += -O3 -g
	NVCCFLAGS  += -O3 -g
endif

DEBUG_LOG ?= 0
ifeq ($(DEBUG_LOG), 1)
	CPPFLAGS   += -DDEBUG_LOG
	NVCCFLAGS  += -DDEBUG_LOG
endif

IEEE ?= 0
ifeq ($(IEEE), 1)
	CPPFLAGS   += -fp-model strict -fp-model source -prec-div -prec-sqrt
	NVCCFLAGS  += -fmad=false -ftz=false -prec-div=true -prec-sqrt=true
endif

CUDA ?= 1
ifeq ($(CUDA), 1)
	CUDA_INC   = -I$(CUDA_INSTALL_PATH)/include
	CUDA_LIB   = -L$(CUDA_INSTALL_PATH)/lib64 -lcudart -lcuda
endif

ALL_INC    += -I$(OPPIC_PATH)/include
ALL_LIBS   += -L$(OPPIC_PATH)/lib

.PHONY: clean mklib

ALL_TARGETS = clean seq omp cuda mpi

all: $(ALL_TARGETS)

mklib:
	@mkdir -p $(OBJ) $(BIN)

common: mklib 

seq: common
		 $(CPP) $(CPPFLAGS) -c advec_op.cpp -o $(OBJ)/advec_op.o $(ALL_INC)
		 $(CPP) $(CPPFLAGS) -c seq/seqkernels.cpp -o $(OBJ)/seqkernels.o $(ALL_INC)
	     $(CPP) $(CPPFLAGS) -o $(BIN)/seq $(OBJ)/advec_op.o $(OBJ)/seqkernels.o $(ALL_INC) $(ALL_LIBS) -loppic_seq

omp: common
	     $(CPP) $(CPPFLAGS) -fopenmp -DUSE_OMP -c advec_op.cpp -o $(OBJ)/advec_op.o $(ALL_INC)
	     $(CPP) $(CPPFLAGS) -fopenmp -DUSE_OMP -c omp/ompkernels.cpp -o $(OBJ)/ompkernels.o $(ALL_INC)
	     $(CPP) $(CPPFLAGS) -fopenmp -o $(BIN)/omp $(OBJ)/advec_op.o $(OBJ)/ompkernels.o $(ALL_INC) $(ALL_LIBS) -loppic_omp

mpi: common
		 $(MPICPP) $(CPPFLAGS) -DUSE_MPI -c advec_op.cpp -o $(OBJ)/advec_op.o $(ALL_INC)
		 $(MPICPP) $(CPPFLAGS) -DUSE_MPI -c mpi/mpikernels.cpp -o $(OBJ)/mpikernels.o $(ALL_INC)
	     $(MPICPP) $(CPPFLAGS) -DUSE_MPI -o $(BIN)/mpi $(OBJ)/advec_op.o $(OBJ)/mpikernels.o $(ALL_INC) $(ALL_LIBS) -loppic_mpi

mpi_hdf5: common
		 $(MPICPP) $(CPPFLAGS) -DUSE_MPI -c advec_op_hdf5.cpp -o $(OBJ)/advec_op_hdf5.o $(ALL_INC)
		 $(MPICPP) $(CPPFLAGS) -DUSE_MPI -c mpi/mpikernels.cpp -o $(OBJ)/mpikernels.o $(ALL_INC)
		 $(MPICPP) $(CPPFLAGS) -DUSE_MPI -c opp_mpi_hdf5.cpp -o $(OBJ)/opp_mpi_hdf5.o $(ALL_INC)
	     $(MPICPP) $(CPPFLAGS) -DUSE_MPI -o $(BIN)/mpi_hdf5 $(OBJ)/advec_op_hdf5.o $(OBJ)/mpikernels.o $(OBJ)/opp_mpi_hdf5.o $(ALL_INC) $(ALL_LIBS) -loppic_mpi -lhdf5


cuda: common
		 $(NVCC) $(NVCCFLAGS) -c cuda/cudakernels.cu -o $(OBJ)/cudakernels.o $(ALL_INC) $(CUDA_INC)
	     $(CPP) $(CPPFLAGS) -o $(BIN)/cuda $(OBJ)/fempic_op.o $(OBJ)/meshes.o $(OBJ)/particles.o $(OBJ)/maths.o $(OBJ)/cudakernels.o $(OBJ)/FESolver.o $(ALL_INC) $(CUDA_INC) $(ALL_LIBS) $(CUDA_LIB) -lopp_cuda

cuda_mpi: mpi_common
		 $(MPICPP) $(CPPFLAGS) -DUSE_MPI -c FESolver.cpp -o $(OBJ)/FESolver.o $(ALL_INC)
		 $(NVCC) $(NVCCFLAGS) -DUSE_MPI -c cuda/cudakernels.cu -o $(OBJ)/cudakernels.o $(ALL_INC) $(CUDA_INC)
	     $(MPICPP) $(CPPFLAGS) -DUSE_MPI -o $(BIN)/cuda_mpi $(OBJ)/fempic_op.o $(OBJ)/meshes.o $(OBJ)/particles.o $(OBJ)/maths.o $(OBJ)/cudakernels.o $(OBJ)/FESolver.o $(ALL_INC) $(CUDA_INC) $(ALL_LIBS) $(CUDA_LIB) -lopp_cuda_mpi

cuda_mpi_hdf5: mpi_common
		 $(NVCC) $(NVCCFLAGS) -DUSE_MPI -c FESolver.cu -o $(OBJ)/FESolver_cu.o $(ALL_INC)	
		 $(MPICPP) $(CPPFLAGS) -DUSE_MPI -c fempic_op_hdf5.cpp -o $(OBJ)/fempic_op_hdf5.o $(ALL_INC)
		 $(MPICPP) $(CPPFLAGS) -DUSE_MPI -c opp_mpi_hdf5.cpp -o $(OBJ)/opp_mpi_hdf5.o $(ALL_INC)
		 $(NVCC) $(NVCCFLAGS) -DUSE_MPI -c cuda/cudakernels.cu -o $(OBJ)/cudakernels.o $(ALL_INC) $(CUDA_INC)
	     $(MPICPP) $(CPPFLAGS) -DUSE_MPI -o $(BIN)/cuda_mpi_hdf5 $(OBJ)/fempic_op_hdf5.o $(OBJ)/opp_mpi_hdf5.o $(OBJ)/meshes.o $(OBJ)/particles.o $(OBJ)/maths.o $(OBJ)/cudakernels.o $(OBJ)/FESolver_cu.o $(ALL_INC) $(CUDA_INC) $(ALL_LIBS) $(CUDA_LIB) -lopp_cuda_mpi -lhdf5

clean: 
	     rm -f *.o *.d *.a
	     rm -f $(OBJ)/*
	     rm -f $(BIN)/*
	     rm -f $(LIB)/*
		 
