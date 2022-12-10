#! /bin/bash

rm -rf files;
mkdir files;

rm -rf fempic

export OMP_NUM_THREADS=4
export OMP_PROC_BIND=close

# Compile Seq without Petsc
mpicxx -o fempic -std=c++11 -Wall -Ofast -fopenmp fempic_op.cpp funcs_from_fempic.cpp field_solver_petsc.cpp seq/seqkernels.cpp ../lib_oppic/oppic.cpp ../lib_oppic/trace.cpp -fPIC -I/opt/homebrew/Cellar/petsc/3.18.1/include -L/opt/homebrew/Cellar/petsc/3.18.1/lib -lpetsc

# Compile Seq with Petsc
# mpicxx -o fempic -std=c++11 -Wall -Ofast -fopenmp fempic_op.cpp funcs_from_fempic.cpp field_solver_original.cpp seq/seqkernels.cpp ../lib_oppic/oppic.cpp ../lib_oppic/trace.cpp -fPIC -I/opt/homebrew/Cellar/petsc/3.18.1/include -L/opt/homebrew/Cellar/petsc/3.18.1/lib -lpetsc


# Compile OMP without Petsc
# mpicxx -o fempic -std=c++17 -Wall -Ofast -fopenmp fempic_op.cpp funcs_from_fempic.cpp field_solver_original.cpp omp/ompkernels.cpp ../lib_oppic/oppic.cpp ../lib_oppic/trace.cpp -fPIC -I/opt/homebrew/Cellar/petsc/3.18.1/include -L/opt/homebrew/Cellar/petsc/3.18.1/lib -lpetsc 

# Compile OMP with Petsc
# mpicxx -o fempic -std=c++17 -Wall -Ofast -fopenmp fempic_op.cpp funcs_from_fempic.cpp field_solver_petsc.cpp omp/ompkernels.cpp ../lib_oppic/oppic.cpp ../lib_oppic/trace.cpp -fPIC -I/opt/homebrew/Cellar/petsc/3.18.1/include -L/opt/homebrew/Cellar/petsc/3.18.1/lib -lpetsc 


# For fempic.cpp prior translation
# mpicxx -o fempic -std=c++17 -Wall -Ofast -fopenmp fempic.cpp funcs_from_fempic.cpp field_solver_petsc.cpp ../lib_oppic/oppic.cpp ../lib_oppic/trace.cpp -fPIC -I/opt/homebrew/Cellar/petsc/3.18.1/include -L/opt/homebrew/Cellar/petsc/3.18.1/lib -lpetsc 

# run fempic_oppic
./fempic
