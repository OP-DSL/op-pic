
# module load gnu-6.3.0/gnu-6.3.0
# module load gnu-6.3.0/ompi-1.10.6
# module load cuda/toolkit-10.2.89

# module load intel-oneapi-2021.3/intel-classic-2021.3 intel-oneapi-2021.3/impi-classic-2021.3
# module load cuda/toolkit-10.2.89

# export LD_LIBRARY_PATH=/ext-home/zl/phd/OP-PIC/app_fempic_oppic/lib:/ext-home/zl/lib_install/petsc-3.18.2-oneapi-2021.3/lib:$LD_LIBRARY_PATH
# export LD_LIBRARY_PATH=/ext-home/zl/phd/OP-PIC/app_fempic_oppic/lib:/ext-home/zl/lib_install/petsc-3.18.2-oneapi-2021.3-release/lib:$LD_LIBRARY_PATH

SRC = .
INC = .
LIB = lib
OBJ = obj
BIN = bin

CPP        = mpicxx
MPICPP	   = mpicxx

CPPFLAGS   = -std=c++11 -Wall -O3 -fopenmp -fPIC
NVCCFLAGS  = -gencode arch=compute_70,code=sm_70 -O3 -m64 -Xptxas -dlcm=ca -Xptxas=-v

DEBUG ?= 0
ifeq ($(DEBUG), 1)
	CPPFLAGS   = -std=c++11 -Wall -O0 -g -fopenmp -fPIC
	NVCCFLAGS  = -gencode arch=compute_70,code=sm_70 -O0 -g -m64 -Xptxas -dlcm=ca -Xptxas=-v
endif

CUDA_INSTALL_PATH = /opt/cuda/10.2.89/toolkit
CUDA_INC   = -I$(CUDA_INSTALL_PATH)/include
CUDA_LIB   = -L$(CUDA_INSTALL_PATH)/lib64 -lcudart -lcuda

# PETSC_PATH = /ext-home/zl/lib_install/petsc-3.18.2-oneapi-2021.3-release
PETSC_PATH = /ext-home/zl/lib_install/petsc-3.18.3-oneapi-2021.3-cuda-rel
PETSC_INC  = -I$(PETSC_PATH)/include
PETSC_LIB  = -L$(PETSC_PATH)/lib -lpetsc

OP2_PATH   = /ext-home/zl/phd/OP2-Common/op2
OP2_INC    = # -I$(OP2_PATH)/include
OP2_LIB    = # -L$(OP2_PATH)/lib

OPPIC_PATH = /ext-home/zl/phd/OP-PIC/lib_oppic
ALL_INC    = -I$(OPPIC_PATH)/include $(PETSC_INC) $(CUDA_INC)
ALL_LIBS   = -L$(OPPIC_PATH)/lib $(PETSC_LIB) $(CUDA_LIB)

.PHONY: clean mklib

# ALL_TARGETS = clean fempic_seq fempic_genseq fempic_openmp fempic_genseq_wopet fempic_openmp_wopet fempic_cuda
ALL_TARGETS = clean test

all: $(ALL_TARGETS)

mklib:
	@mkdir -p $(OBJ) $(BIN)

test: 
	     $(CPP) $(CPPFLAGS) -c fempic_op.cpp -o $(OBJ)/fempic_op.o $(ALL_INC)
	     $(CPP) $(CPPFLAGS) -c fempic_ori/meshes.cpp -o $(OBJ)/meshes.o $(ALL_INC)
	     $(CPP) $(CPPFLAGS) -c fempic_ori/particles.cpp -o $(OBJ)/particles.o $(ALL_INC)
		 $(CPP) $(CPPFLAGS) -c fempic_ori/maths.cpp -o $(OBJ)/maths.o $(ALL_INC)
		 $(CPP) $(CPPFLAGS) -c FESolver.cpp -o $(OBJ)/FESolver.o $(ALL_INC)

seq: mklib test
		 $(CPP) $(CPPFLAGS) -c seq/seqkernels.cpp -o $(OBJ)/seqkernels.o $(ALL_INC)
	     $(CPP) $(CPPFLAGS) -o $(BIN)/seq $(OBJ)/fempic_op.o $(OBJ)/meshes.o $(OBJ)/particles.o $(OBJ)/maths.o $(OBJ)/seqkernels.o $(OBJ)/FESolver.o $(ALL_INC) $(ALL_LIBS) -loppic_seq

omp: mklib test
		 $(CPP) $(CPPFLAGS) -c omp/ompkernels.cpp -o $(OBJ)/ompkernels.o $(ALL_INC)
	     $(CPP) $(CPPFLAGS) -o $(BIN)/omp $(OBJ)/fempic_op.o $(OBJ)/meshes.o $(OBJ)/particles.o $(OBJ)/maths.o $(OBJ)/ompkernels.o $(OBJ)/FESolver.o $(ALL_INC) $(ALL_LIBS) -loppic_seq


clean: 
	     rm -f *.o *.d *.a
	     rm -f $(OBJ)/*
	     rm -f $(BIN)/*
	     rm -f $(LIB)/*
		 
# fempic_seq: mklib
# 	     $(CPP) $(CPPFLAGS) -c fempic.cpp -o $(OBJ)/fempic.o $(ALL_INC)
# 	     $(CPP) $(CPPFLAGS) -c funcs_from_fempic.cpp -o $(OBJ)/funcs_from_fempic.o $(ALL_INC)
# 	     $(CPP) $(CPPFLAGS) -c field_solver_original.cpp -o $(OBJ)/field_solver_original.o $(ALL_INC)
# 	     $(CPP) $(CPPFLAGS) -o $(BIN)/fempic_seq $(OBJ)/fempic.o $(OBJ)/funcs_from_fempic.o $(OBJ)/field_solver_original.o $(ALL_INC) $(ALL_LIBS) -loppic_seq

# fempic_genseq: mklib
# 	     $(CPP) $(CPPFLAGS) -c fempic_op.cpp -o $(OBJ)/fempic_op.o $(ALL_INC)
# 	     $(CPP) $(CPPFLAGS) -c funcs_from_fempic.cpp -o $(OBJ)/funcs_from_fempic.o $(ALL_INC)
# 	     $(CPP) $(CPPFLAGS) -c field_solver_petsc.cpp -o $(OBJ)/field_solver_petsc.o $(ALL_INC)
# 	     $(CPP) $(CPPFLAGS) -c seq/seqkernels.cpp -o $(OBJ)/seqkernels.o $(ALL_INC)
# 	     $(CPP) $(CPPFLAGS) -o $(BIN)/fempic_genseq $(OBJ)/fempic_op.o $(OBJ)/funcs_from_fempic.o $(OBJ)/field_solver_petsc.o $(OBJ)/seqkernels.o $(ALL_INC) $(ALL_LIBS) -loppic_seq 

# fempic_openmp: mklib
# 	     $(CPP) $(CPPFLAGS) -c fempic_op.cpp -o $(OBJ)/fempic_op.o $(ALL_INC)
# 	     $(CPP) $(CPPFLAGS) -c funcs_from_fempic.cpp -o $(OBJ)/funcs_from_fempic.o $(ALL_INC)
# 	     $(CPP) $(CPPFLAGS) -c field_solver_petsc.cpp -o $(OBJ)/field_solver_petsc.o $(ALL_INC)
# 	     $(CPP) $(CPPFLAGS) -c omp/ompkernels.cpp -o $(OBJ)/ompkernels.o $(ALL_INC)
# 	     $(CPP) $(CPPFLAGS) -o $(BIN)/fempic_openmp $(OBJ)/fempic_op.o $(OBJ)/funcs_from_fempic.o $(OBJ)/field_solver_petsc.o $(OBJ)/ompkernels.o $(ALL_INC) $(ALL_LIBS) -loppic_omp 

# fempic_genseq_wopet: mklib
# 	     $(CPP) $(CPPFLAGS) -c fempic_op.cpp -o $(OBJ)/fempic_op.o $(ALL_INC)
# 	     $(CPP) $(CPPFLAGS) -c funcs_from_fempic.cpp -o $(OBJ)/funcs_from_fempic.o $(ALL_INC)
# 	     $(CPP) $(CPPFLAGS) -c field_solver_original.cpp -o $(OBJ)/field_solver_original.o $(ALL_INC)
# 	     $(CPP) $(CPPFLAGS) -c seq/seqkernels.cpp -o $(OBJ)/seqkernels.o $(ALL_INC)
# 	     $(CPP) $(CPPFLAGS) -o $(BIN)/fempic_genseq_wopet $(OBJ)/fempic_op.o $(OBJ)/funcs_from_fempic.o $(OBJ)/field_solver_original.o $(OBJ)/seqkernels.o  $(ALL_INC) $(ALL_LIBS) -loppic_seq 

# fempic_openmp_wopet: mklib
# 	     $(CPP) $(CPPFLAGS) -c fempic_op.cpp -o $(OBJ)/fempic_op.o $(ALL_INC)
# 	     $(CPP) $(CPPFLAGS) -c funcs_from_fempic.cpp -o $(OBJ)/funcs_from_fempic.o $(ALL_INC)
# 	     $(CPP) $(CPPFLAGS) -c field_solver_original.cpp -o $(OBJ)/field_solver_original.o $(ALL_INC)
# 	     $(CPP) $(CPPFLAGS) -c omp/ompkernels.cpp -o $(OBJ)/ompkernels.o $(ALL_INC)
# 	     $(CPP) $(CPPFLAGS) -o $(BIN)/fempic_openmp_wopet $(OBJ)/fempic_op.o $(OBJ)/funcs_from_fempic.o $(OBJ)/field_solver_original.o $(OBJ)/ompkernels.o $(ALL_INC) $(ALL_LIBS) -loppic_omp 

# fempic_cuda: mklib
# 	     $(CPP) $(CPPFLAGS) -c fempic_op.cpp -o $(OBJ)/fempic_op.o $(ALL_INC)
# 	     $(CPP) $(CPPFLAGS) -c funcs_from_fempic.cpp -o $(OBJ)/funcs_from_fempic.o $(ALL_INC)
# 	     $(CPP) $(CPPFLAGS) -c cuda/field_solver_petsc_cuda.cpp -o $(OBJ)/field_solver_petsc_cuda.o $(ALL_INC)
# 	     nvcc $(NVCCFLAGS) -c cuda/cudakernels.cu -o $(OBJ)/cudakernels.o $(ALL_INC)
# 		 $(CPP) $(CPPFLAGS) -o $(BIN)/fempic_cuda $(OBJ)/fempic_op.o $(OBJ)/funcs_from_fempic.o $(OBJ)/field_solver_petsc_cuda.o $(OBJ)/cudakernels.o $(ALL_INC) $(ALL_LIBS) -loppic_cuda -lcudart

# fempic_cuda_wopet: mklib
# 	     $(CPP) $(CPPFLAGS) -c fempic_op.cpp -o $(OBJ)/fempic_op.o $(ALL_INC)
# 	     $(CPP) $(CPPFLAGS) -c funcs_from_fempic.cpp -o $(OBJ)/funcs_from_fempic.o $(ALL_INC)
# 	     $(CPP) $(CPPFLAGS) -c cuda/field_solver_original_cuda.cpp -o $(OBJ)/field_solver_original_cuda.o $(ALL_INC)
# 	     nvcc $(NVCCFLAGS) -c cuda/cudakernels.cu -o $(OBJ)/cudakernels.o $(ALL_INC)
# 		 $(CPP) $(CPPFLAGS) -o $(BIN)/fempic_cuda $(OBJ)/fempic_op.o $(OBJ)/funcs_from_fempic.o $(OBJ)/field_solver_original_cuda.o $(OBJ)/cudakernels.o $(ALL_INC) $(ALL_LIBS) -loppic_cuda -lcudart



