
# module load gnu-6.3.0/gnu-6.3.0
# module load gnu-6.3.0/ompi-1.10.6
# module load cuda/toolkit-10.2.89

CPP        = mpicxx
CPPFLAGS   = -std=c++11 -Wall -O3 -fopenmp -fPIC
MPICPP	   = mpicxx
MPIFLAGS   = -std=c++11 -Wall -O3 -fopenmp -fPIC
NVCCFLAGS  = -O3 -m64 -Xptxas -dlcm=ca -Xptxas=-v -use_fast_math -arch=compute_60 -code=sm_60

CUDA_INSTALL_PATH = # /opt/cuda/10.2.89/toolkit
CUDA_INC   = -I$(CUDA_INSTALL_PATH)/include
CUDA_LIB   = -L$(CUDA_INSTALL_PATH)/lib64

PETSC_INC  = # -I/opt/homebrew/Cellar/petsc/3.18.1/include
PETSC_LIB  = # -L/opt/homebrew/Cellar/petsc/3.18.1/lib -lpetsc

ALL_TARGETS = clean_all oppic_seq fempic_seq fempic_genseq fempic_openmp fempic_genseq_wopet fempic_openmp_wopet fempic_cuda

all: $(ALL_TARGETS)

oppic_seq: oppic.o trace.o
	     $(CPP) $(CPPFLAGS) -shared -I../lib_oppic -o liboppic_seq.so $^

oppic.o: ../lib_oppic/oppic.cpp
	     $(CPP) $(CPPFLAGS) -c $^ -o $@

trace.o: ../lib_oppic/trace.cpp
	     $(CPP) $(CPPFLAGS) -c $^ -o $@



fempic_seq: fempic.o funcs_from_fempic.o field_solver_petsc.o
	     $(CPP) $(CPPFLAGS) $(PETSC_INC) $(PETSC_LIB) -L. -loppic_seq $^ -o $@

fempic_genseq: fempic_op.o funcs_from_fempic.o field_solver_petsc.o seq/seqkernels.o
	     $(CPP) $(CPPFLAGS) $(PETSC_INC) $(PETSC_LIB) -L. -loppic_seq $^ -o $@

fempic_openmp: fempic_op.o funcs_from_fempic.o field_solver_petsc.o omp/ompkernels.o
	     $(CPP) $(CPPFLAGS) $(PETSC_INC) $(PETSC_LIB) -L. -loppic_seq $^ -o $@

fempic_genseq_wopet: fempic_op.o funcs_from_fempic.o field_solver_original.o seq/seqkernels.o
	     $(CPP) $(CPPFLAGS) $(PETSC_INC) $(PETSC_LIB) -L. -loppic_seq $^ -o $@

fempic_openmp_wopet: fempic_op.o funcs_from_fempic.o field_solver_original.o omp/ompkernels.o
	     $(CPP) $(CPPFLAGS) $(PETSC_INC) $(PETSC_LIB) -L. -loppic_seq $^ -o $@

fempic_cuda: fempic_op.o funcs_from_fempic.o field_solver_original.o cudakernels.o
	     $(CPP) $(CPPFLAGS) $(PETSC_INC) $(PETSC_LIB) $(CUDA_INC) $(CUDA_LIB) -I. -L. -loppic_seq $^ -lcudart -o $@




%:
	     $(CPP) $(CPPFLAGS) $(PETSC_INC) $(PETSC_LIB) -c %.cpp -o %

cudakernels.o:	
		nvcc $(VAR) $(INC) $(NVCCFLAGS) -c cuda/cudakernels.cu



clean_all:
		rm -f fempic_seq fempic_genseq fempic_openmp fempic_genseq_wopet fempic_openmp_wopet fempic_cuda *.o cuda/*.o *.so

clean: 
		rm -f fempic_seq fempic_genseq fempic_openmp fempic_genseq_wopet fempic_openmp_wopet fempic_cuda *.o cuda/*.o

