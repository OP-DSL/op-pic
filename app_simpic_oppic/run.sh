#!/bin/bash

rm -rf files;
mkdir files;

rm -rf fempic

export OMP_NUM_THREADS=4
export OMP_PROC_BIND=close

# module load gnu-9.3.0/openmpi-4.0.4
# export LD_LIBRARY_PATH=/ext-home/zl/phd/OP-PIC/app_simpic_oppic:$LD_LIBRARY_PATH

NT=50               # of time steps
PPC=100             # of particles per cell
CPP=10              # of cell per process
DTFACTOR=0.001      #defines a fraction of dt vs. time steps from plasma frequency;
LHSV=20000          #applied voltage on left-hand side; RHS is grounded; 

CMDLINE3="./simpic_seq -ppc $PPC -ncpp $CPP -nt $NT -dtfactor $DTFACTOR -lhsv $LHSV"
echo $CMDLINE3
$CMDLINE3

CMDLINE3="./simpic_genseq -ppc $PPC -ncpp $CPP -nt $NT -dtfactor $DTFACTOR -lhsv $LHSV"
echo $CMDLINE3
$CMDLINE3

CMDLINE3="./simpic_openmp -ppc $PPC -ncpp $CPP -nt $NT -dtfactor $DTFACTOR -lhsv $LHSV"
echo $CMDLINE3
$CMDLINE3

echo "Completed..."

