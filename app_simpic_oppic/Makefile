
CPP        = mpicxx
CPPFLAGS   = -std=c++11 -Wall -O3 -fopenmp -fPIC
MPICPP	   = mpicxx
MPIFLAGS   = -std=c++11 -Wall -O3 -fopenmp -fPIC

PETSC_INC  = # -I/opt/homebrew/Cellar/petsc/3.18.1/include
PETSC_LIB  = # -L/opt/homebrew/Cellar/petsc/3.18.1/lib -lpetsc

ALL_TARGETS = clean_all oppic_seq simpic_seq simpic_genseq simpic_openmp

all: $(ALL_TARGETS)

oppic_seq: oppic.o trace.o
	     $(CPP) $(CPPFLAGS) -shared -I../lib_oppic -o liboppic_seq.so $^

oppic.o: ../lib_oppic/oppic.cpp
	     $(CPP) $(CPPFLAGS) -c $^ -o $@

trace.o: ../lib_oppic/trace.cpp
	     $(CPP) $(CPPFLAGS) -c $^ -o $@



simpic_seq: simpic.cpp funcs_from_simpic.cpp
	     $(CPP) $(CPPFLAGS) $(PETSC_INC) $(PETSC_LIB) -L. -loppic_seq $^ -o $@

simpic_genseq: simpic_op.cpp funcs_from_simpic.cpp seq/seqkernels.cpp
	     $(CPP) $(CPPFLAGS) $(PETSC_INC) $(PETSC_LIB) -L. -loppic_seq $^ -o $@

simpic_openmp: simpic_op.cpp funcs_from_simpic.cpp omp/ompkernels.cpp
	     $(CPP) $(CPPFLAGS) $(PETSC_INC) $(PETSC_LIB) -L. -loppic_seq $^ -o $@

simpic_cuda: simpic_op.cpp funcs_from_simpic.cpp seq/seqkernels.cpp
	     $(CPP) $(CPPFLAGS) $(PETSC_INC) $(PETSC_LIB) -L. -loppic_seq $^ -o $@

clean_all:
		rm -f simpic_seq simpic_genseq simpic_openmp *.o cuda/*.o *.so

clean: 
		rm -f simpic_seq simpic_genseq simpic_openmp *.o cuda/*.o

