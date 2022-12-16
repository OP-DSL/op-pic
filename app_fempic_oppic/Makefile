
CPP        = mpicxx
CPPFLAGS   = -std=c++11 -Wall -O3 -fopenmp -fPIC
MPICPP	   = mpicxx
MPIFLAGS   = -std=c++11 -Wall -O3 -fopenmp -fPIC

PETSC_INC  = -I/opt/homebrew/Cellar/petsc/3.18.1/include
PETSC_LIB  = -L/opt/homebrew/Cellar/petsc/3.18.1/lib -lpetsc

ALL_TARGETS = clean_all oppic_seq fempic_seq fempic_genseq fempic_openmp fempic_genseq_wopet fempic_openmp_wopet

all: $(ALL_TARGETS)

oppic_seq: oppic.o trace.o
	     $(CPP) $(CPPFLAGS) -shared -I../lib_oppic -o liboppic_seq.so $^

oppic.o: ../lib_oppic/oppic.cpp
	     $(CPP) $(CPPFLAGS) -c $^ -o $@

trace.o: ../lib_oppic/trace.cpp
	     $(CPP) $(CPPFLAGS) -c $^ -o $@



fempic_seq: fempic.cpp funcs_from_fempic.cpp field_solver_petsc.cpp
	     $(CPP) $(CPPFLAGS) $(PETSC_INC) $(PETSC_LIB) -L. -loppic_seq $^ -o $@

fempic_genseq: fempic_op.cpp funcs_from_fempic.cpp field_solver_petsc.cpp seq/seqkernels.cpp
	     $(CPP) $(CPPFLAGS) $(PETSC_INC) $(PETSC_LIB) -L. -loppic_seq $^ -o $@

fempic_openmp: fempic_op.cpp funcs_from_fempic.cpp field_solver_petsc.cpp omp/ompkernels.cpp
	     $(CPP) $(CPPFLAGS) $(PETSC_INC) $(PETSC_LIB) -L. -loppic_seq $^ -o $@

fempic_genseq_wopet: fempic_op.cpp funcs_from_fempic.cpp field_solver_original.cpp seq/seqkernels.cpp
	     $(CPP) $(CPPFLAGS) $(PETSC_INC) $(PETSC_LIB) -L. -loppic_seq $^ -o $@

fempic_openmp_wopet: fempic_op.cpp funcs_from_fempic.cpp field_solver_original.cpp omp/ompkernels.cpp
	     $(CPP) $(CPPFLAGS) $(PETSC_INC) $(PETSC_LIB) -L. -loppic_seq $^ -o $@

fempic_cuda: fempic_op.cpp funcs_from_fempic.cpp field_solver_petsc.cpp seq/seqkernels.cpp
	     $(CPP) $(CPPFLAGS) $(PETSC_INC) $(PETSC_LIB) -L. -loppic_seq $^ -o $@

clean_all:
		rm -f fempic_seq fempic_genseq fempic_openmp fempic_genseq_wopet fempic_openmp_wopet *.o cuda/*.o *.so

clean: 
		rm -f fempic_seq fempic_genseq fempic_openmp fempic_genseq_wopet fempic_openmp_wopet *.o cuda/*.o

