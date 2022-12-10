#!/bin/bash

rm -f simpic files/*.dat

# Compile Seq
mpicxx -o simpic -std=c++11 -Wall -Ofast -fopenmp simpic_op.cpp funcs_from_simpic.cpp seq/seqkernels.cpp ../lib_oppic/oppic.cpp ../lib_oppic/trace.cpp 

# Compile OMP
# mpicxx -o simpic -std=c++11 -Wall -Ofast -fopenmp simpic_op.cpp funcs_from_simpic.cpp omp/ompkernels.cpp ../lib_oppic/oppic.cpp ../lib_oppic/trace.cpp 

# For simpic.cpp prior translation
# mpicxx -o simpic -std=c++11 -Wall -Ofast -fopenmp simpic.cpp funcs_from_simpic.cpp ../lib_oppic/oppic.cpp ../lib_oppic/trace.cpp 


NT=50               # of time steps
PPC=100             # of particles per cell
CPP=10              # of cell per process
DTFACTOR=0.001      #defines a fraction of dt vs. time steps from plasma frequency;
LHSV=20000          #applied voltage on left-hand side; RHS is grounded; 

CMDLINE3="./simpic -ppc $PPC -ncpp $CPP -nt $NT -dtfactor $DTFACTOR -lhsv $LHSV"
echo $CMDLINE3
$CMDLINE3

echo "Completed..."

